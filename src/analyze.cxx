#include "analyze.hxx"
#include "analyze/bond_orders.hxx"
#include "append.hxx"
#include "elements.hxx"
#include "graph.hxx"
#include "geom.hxx"
#include "clone.hxx"
#include "contacts.hxx"
#include "pfx/pfx.hxx"
#include "smiles.hxx"
#include <numeric>
#include <queue>
#include <stdio.h>
#include <tuple>

#if defined(WIN32) && !defined(drand48)
#define drand48() ((double)rand()/(double)RAND_MAX)
#endif
#if defined(WIN32) && !defined(srand48)
#define srand48 srand
#endif

namespace {
    using namespace desres::msys;
    struct BondFinder {
        SystemPtr mol;
        BondFinder(SystemPtr m) : mol(m) {}

        /* Excludes H-H and pseudo-pseudo bonds */
        bool exclude(Id i, Id j) const {
            if (i>=j) return true;
            int ai = mol->atomFAST(i).atomic_number;
            int aj = mol->atomFAST(j).atomic_number;
            if ((ai==1 && aj==1) || (ai==0 && aj==0))
                return true;
            return false;
        }

        void operator()(Id i, Id j, double d2) const {
            if (d2>0.00001) { // some virtuals can be very close to the host
                int ai = mol->atomFAST(i).atomic_number;
                int aj = mol->atomFAST(j).atomic_number;
                double ri = RadiusForElement(ai);
                double rj = RadiusForElement(aj);
                double cut = 0.6 * (ri+rj);
                if (d2 < cut*cut) {
                    mol->addBond(i,j);
                }
            }
        }
    };
}

namespace desres { namespace msys {

    SystemPtr CanonicalizeMoleculeByTopids(SystemPtr mol, IdList const& atoms,
                                           std::map<Id, Id> &aid_to_canId,
                                           std::map<Id, Id> &bid_to_canId){
        aid_to_canId.clear();
        bid_to_canId.clear();

        IdList topids = ComputeTopologicalIds(mol);
        std::map<Id, IdList> topid_to_atomid;
        for(Id aid : atoms){
            topid_to_atomid[topids.at(aid)].push_back(aid);
        }
        IdList neworder;
        for (auto const& kv : topid_to_atomid) {
            for( auto const& v: kv.second){
                aid_to_canId[v] = neworder.size();
                neworder.push_back(v);
            }
        }
        std::vector< std::tuple<Id, Id, Id, Id, Id> > newBonds;
        for (Id bid : mol->bonds()) {
            bond_t const& bond = mol->bond(bid);
            auto it0 = aid_to_canId.find(bond.i);
            auto it1 = aid_to_canId.find(bond.j);
            /* only keep the bond if we have both atoms */
            if(it0 == aid_to_canId.end() or it1 == aid_to_canId.end()) continue;
            Id newAtom0 = it0->second;
            Id newAtom1 = it1->second;
            Id tid0 = topids.at(bond.i);
            Id tid1 = topids.at(bond.j);
            if ( tid1 < tid0 ){
                std::swap(tid0, tid1);
            }
            if ( newAtom1 < newAtom0){
                std::swap(newAtom0, newAtom1);
            }
            newBonds.push_back(std::make_tuple(tid0, tid1, newAtom0, newAtom1, bid));
        }
        std::sort(newBonds.begin(), newBonds.end());
        SystemPtr newmol = Clone(mol, neworder);
        /* remove all the old bonds */
            for (Id bid : newmol->bonds()){
                newmol->delBond(bid);
        }
        /* we can clone again here if we want the bonds to start back at id=0 */
        // newmol = Clone(newmol, neworder);

        /* add in the new bonds and keep track of their ids */
        for (auto const& v : newBonds){
            Id canid = newmol->addBond(std::get<2>(v), std::get<3>(v));
            bid_to_canId[std::get<4>(v)] = canid;
        }
        return newmol;
    }


    void AssignBondOrderAndFormalCharge(SystemPtr mol, unsigned flags, std::chrono::milliseconds timeout) {
        auto deadline = std::chrono::system_clock::now() + timeout;
        MultiIdList fragments;
        mol->updateFragids(&fragments);
        IdList pmap(mol->maxAtomId(), BadId);

        /* will compute graphs lazily */
        std::vector<GraphPtr> graphs(fragments.size());
        typedef std::map<std::string, IdList> FragmentHash;
        FragmentHash fragment_hash;

        for (Id i=0; i<fragments.size(); i++) {
            fragment_hash[Graph::hash(mol, fragments[i])].push_back(i);
        }
        FragmentHash::iterator it;
        for (it=fragment_hash.begin(); it!=fragment_hash.end(); ++it) {
            /* fetch fragments with the same formula */
            IdList& frags = it->second;
            /* unique formula -> unique fragment */
            if (frags.size()==1) {
                AssignBondOrderAndFormalCharge(mol, fragments[frags[0]], INT_MAX, flags, std::chrono::duration_cast<std::chrono::milliseconds>(deadline - std::chrono::system_clock::now()));
                continue;
            }

            /* We have multiple fragments with the same formula.  */
            for (Id frag : frags) {
                graphs[frag] = Graph::create(mol, fragments[frag]);
            }
            std::vector<IdPair> perm;
            while (!frags.empty()) {
                AssignBondOrderAndFormalCharge(mol, fragments[frags[0]], INT_MAX, flags, std::chrono::duration_cast<std::chrono::milliseconds>(deadline - std::chrono::system_clock::now()));
                IdList unmatched;
                GraphPtr ref = graphs[frags[0]];
                for (Id i=1; i<frags.size(); i++) {
                    GraphPtr sel = graphs[frags[i]];
                    if (!ref->match(sel, perm)) {
                        /* didn't match, so push into the next iteration
                         * for another bond order calculation */
                        unmatched.push_back(frags[i]);
                    } else {
                        /* map atom properties */
                        for (IdPair const&p : perm) {
                            const Id ai = p.first;
                            const Id bi = p.second;
                            mol->atom(bi).formal_charge = mol->atom(ai).formal_charge;
                            pmap.at(ai) = bi;
                        }
                        /* map bond properties */
                        for (IdPair const&p : perm) {
                            const Id ai = p.first;
                            const Id bi = p.second;
                            for (Id bnd : mol->bondsForAtom(ai)) {
                                bond_t const& src = mol->bond(bnd);
                                const Id aj = src.other(ai);
                                if (ai>aj) continue;
                                const Id bj = pmap.at(aj);
                                if (bad(bj)) continue;
                                bond_t& dst = mol->bond(mol->findBond(bi,bj));
                                dst.order = src.order;
                            }
                        }
                    }
                }
                frags.swap(unmatched);
            }
        }
    }


    void AssignBondOrderAndFormalCharge(SystemPtr mol,
                                        IdList const& atoms,
                                        int total_charge,
                                        unsigned flags,
                                        std::chrono::milliseconds timeout) {
#ifdef MSYS_WITHOUT_LPSOLVE
        MSYS_FAIL("LPSOLVE functionality was not included.");
#else
        if (atoms.empty()) return;

        /* get a canonicalized ordering for the molecule based on topids */
        std::map<Id, Id> aid_to_canId;
        std::map<Id, Id> bid_to_canId;
        SystemPtr canmol = CanonicalizeMoleculeByTopids(mol, atoms, aid_to_canId, bid_to_canId);

        bool compute_resonant_charge = flags & AssignBondOrder::ComputeResonantCharges;
        BondOrderAssigner boa(canmol, canmol->atoms(), compute_resonant_charge, timeout);
        if (total_charge != INT_MAX) {
            boa.setTotalCharge(total_charge);
        }
        boa.solveIntegerLinearProgram();
        boa.assignSolutionToAtoms();

        Id qprop = BadId, oprop = BadId;
        Id can_qprop = BadId, can_oprop = BadId;
        if (compute_resonant_charge) {
            qprop = mol->addAtomProp("resonant_charge", FloatType);
            oprop = mol->addBondProp("resonant_order", FloatType);
            can_qprop = canmol->atomPropIndex("resonant_charge");
            can_oprop = canmol->bondPropIndex("resonant_order");
        }
        /* copy fc/bo from canonical mol to input mol */
        for ( auto const& kv : aid_to_canId){
            mol->atom(kv.first).formal_charge = canmol->atom(kv.second).formal_charge;
            if (compute_resonant_charge) {
                mol->atomPropValue(kv.first, qprop) = canmol->atomPropValue(kv.second, can_qprop);
            }
        }
        for ( auto const& kv : bid_to_canId){
            mol->bond(kv.first).order = canmol->bond(kv.second).order;
            if (compute_resonant_charge) {
                mol->bondPropValue(kv.first, oprop) = canmol->bondPropValue(kv.second, can_oprop);
            }
        }
#endif
    }


    void GuessBondConnectivity(SystemPtr mol, bool periodic) {
        const IdList atoms(mol->atoms());
        static const double cutoff = 4.0;
        std::vector<Float> pos(3*mol->maxAtomId());
        if (pos.empty()) return;
        for (Id i : atoms) {
            atom_t const& atom = mol->atom(i);
            pos[3*i  ] = atom.x;
            pos[3*i+1] = atom.y;
            pos[3*i+2] = atom.z;
        }

        if (periodic) {
            /* brute force method that handles triclinics */
            const double* cell = mol->global_cell[0];
            Float box[9], inv[9], vec[3];
            pfx::trans_3x3(box, cell);
            if (!pfx::inverse_3x3(inv, box)) {
                memset(inv, 0, sizeof(inv));
            }
            std::copy(cell, cell+9, box);

            for (Id i=0, n=atoms.size(); i<n; i++) {
                const auto& iatm = mol->atomFAST(atoms[i]);
                const auto ai = iatm.atomic_number;
                const auto ri = RadiusForElement(ai);
                for (Id j=i+1; j<n; j++) {
                    const auto& jatm = mol->atomFAST(atoms[j]);
                    const auto aj = jatm.atomic_number;
                    // don't bond H-H or Virt-Virt
                    if ((ai==1 && aj==1) || (ai==0 && aj==0)) continue;
                    const auto rj = RadiusForElement(aj);
                    const double cut = 0.6 * (ri+rj);
                    Float d[3] = {jatm.x-iatm.x, jatm.y-iatm.y, jatm.z-iatm.z};
                    std::copy(d,d+3,vec);
                    pfx::wrap_vector(box, inv, d);
                    Float dx = d[0]+vec[0];
                    Float dy = d[1]+vec[1];
                    Float dz = d[2]+vec[2];
                    Float d2 = dx*dx + dy*dy + dz*dz;
                    if (d2 < cut*cut) {
                        mol->addBond(atoms[i], atoms[j]);
                    }
                }
            }

        } else {
            BondFinder finder(mol);
            if (mol->ctCount()==1) {
                find_contacts(cutoff, &pos[0],
                              atoms.begin(), atoms.end(),
                              atoms.begin(), atoms.end(),
                              finder);
            } else {
                for (Id i=0, n=mol->ctCount(); i<n; i++) {
                    IdList const& atoms = mol->atomsForCt(i);
                    find_contacts(cutoff, &pos[0],
                                  atoms.begin(), atoms.end(),
                                  atoms.begin(), atoms.end(),
                                  finder);
                }
            }
        }
        /* if a pseudo has multiple bonds:
              keep the shortest non-pseudo bond (ie closest atom is "parent"),
           if a hydrogen has multiple bonds:
              keep the shortest heavy bond (preserving any virtual bonds)
        */
        for (Id i=0; i<mol->maxAtomId(); i++) {
            if (!mol->hasAtom(i)) continue;
            int anum=mol->atomFAST(i).atomic_number;
            if (anum>1) continue;
            if (mol->bondCountForAtom(i)<=1) continue;
            IdList candidates;
            for (Id b : mol->bondsForAtom(i)){
                int anum2=mol->atomFAST(mol->bond(b).other(i)).atomic_number;
                if(anum2 <= anum) continue;
                candidates.push_back(b);
            }
            if(candidates.size()<=1)continue;
            Id shortest_bond = BadId;
            double shortest_dist = HUGE_VAL;
            const double* pi = &pos[3*i];
            const double x=pi[0];
            const double y=pi[1];
            const double z=pi[2];
            for (Id b : candidates) {
                Id j = mol->bond(b).other(i);
                const double* pj = &pos[3*j];
                const double dx = pj[0]-x;
                const double dy = pj[1]-y;
                const double dz = pj[2]-z;
                const double d2 = dx*dx + dy*dy + dz*dz;
                if (d2<shortest_dist) {
                    shortest_bond = b;
                    shortest_dist = d2;
                }
            }
            assert(!bad(shortest_bond));
            for (Id b : candidates) {
                if (b!=shortest_bond) mol->delBond(b);
            }
        }
    }

    std::map<Id,IdList> FindDistinctFragments(SystemPtr mol, MultiIdList const& fragments, std::vector<std::string> const& keys) {
        std::map<Id, IdList> result;
        /* will compute graphs lazily */
        std::vector<GraphPtr> graphs(fragments.size());
        typedef std::map<std::string, IdList> FragmentHash;
        FragmentHash fragment_hash;
        for (Id i=0; i<fragments.size(); i++) {
            std::string key;
            if (keys.empty()) {
                key = Graph::hash(mol, fragments[i]);
            } else {
                key = keys.at(i);
            }
            fragment_hash[key].push_back(i);
        }
        FragmentHash::iterator it;
        for (it=fragment_hash.begin(); it!=fragment_hash.end(); ++it) {
            /* unique formula -> unique fragment */
            IdList& frags = it->second;
            if (frags.size()==1) {
                result[frags[0]] = frags;
                continue;
            }
            /* must do isomorphism checks. */
            for (Id frag : frags) {
                graphs[frag] = Graph::create(mol, fragments[frag]);
            }
            std::vector<IdPair> perm;
            while (!frags.empty()) {
                Id fragid = frags[0];
                IdList& matched = result[fragid];
                matched.push_back(fragid);
                IdList unmatched;
                GraphPtr ref = graphs[fragid];
                for (Id i=1; i<frags.size(); i++) {
                    GraphPtr sel = graphs[frags[i]];
                    if (!ref->match(sel, perm)) {
                        unmatched.push_back(frags[i]);
                    } else {
                        matched.push_back(frags[i]);
                    }
                }
                frags.swap(unmatched);
            }
        }
        return result;
    }

}}

namespace {

    bool has_water_residue_name( const std::string& resname ) {
        static const char * names[] = {
            "H2O", "HH0", "OHH", "HOH", "OH2", "SOL", "WAT", "TIP",
            "TIP2", "TIP3", "TIP4", "SPC"
        };
        unsigned i,n = sizeof(names)/sizeof(names[0]);
        for (i=0; i<n; i++) {
            if (resname==names[i]) return true;
        }
        return false;
    }

    void find_sidechain(System* mol, Id res, Id ca) {
        /* pick a c-beta atom, or settle for a hydrogen */
        Id cb=BadId;
        for (Id nbr : mol->bondedAtoms(ca)) {
            if (mol->atomFAST(nbr).type==AtomProBack) continue;
            if (bad(cb)) {
                cb=nbr;
            } else {
                /* already have a cb candidate.  Pick the better one. */
                if (mol->atomFAST(nbr).atomic_number==1 ||
                    toupper(mol->atomFAST(nbr).name[0]=='H')) continue;
                cb=nbr;
            }
        }
        if (bad(cb)) return;
        /* starting from cb, recursively add all atoms in the residue
         * which are bonded to the cb but are not backbone. */
        std::queue<Id> q;
        q.push(cb);
        while (!q.empty()) {
            Id atm = q.front();
            mol->atomFAST(atm).type = AtomProSide;
            for (Id nbr : mol->bondedAtoms(atm)) {
                atom_t const& nbratm = mol->atomFAST(nbr);
                if (nbratm.type==AtomOther && nbratm.residue==res) {
                    q.push(nbr);
                }
            }
            q.pop();
        }
    }

    // FIXME: I've tried writing a water topology check by hand, but it's
    // always been slower than using Graph::match.  Go figure.
    static std::string waterhash;
    static GraphPtr watergraph;
    static std::string ace_hash;
    static GraphPtr ace_graph;
    static std::string nme_hash;
    static GraphPtr nme_graph;
    struct _ {
        _() {
            make_watergraph();
            make_ace_graph();
            make_nme_graph();
        }
        void make_watergraph() {
            SystemPtr m = System::create();
            m->addChain();
            m->addResidue(0);
            m->addAtom(0);
            m->addAtom(0);
            m->addAtom(0);
            m->atom(0).atomic_number=8;
            m->atom(1).atomic_number=1;
            m->atom(2).atomic_number=1;
            m->addBond(0,1);
            m->addBond(0,2);
            watergraph = Graph::create(m,m->atoms());
            waterhash = Graph::hash(m,m->atoms());
        }
        void make_ace_graph() {
            auto m = FromSmilesString("CC(=O)N");
            ace_graph = Graph::create(m, {0,1,2,4,5,6});
            ace_hash = Graph::hash(   m, {0,1,2,4,5,6});
        }
        void make_nme_graph() {
            auto m = FromSmilesString("CN");
            nme_graph = Graph::create(m, {0,1,2,3,4,5});
            nme_hash = Graph::hash(   m, {0,1,2,3,4,5});
        }
    } static_initializer;

    typedef std::map<std::string,AtomType> NameMap;

    static bool is_protein_capping_group(SystemPtr mol, IdList const& atoms) {
        std::vector<IdPair> perm;
        GraphPtr g;
        auto h = Graph::hash(mol, atoms);
        if (h==ace_hash) {
            if (!g) g = Graph::create(mol, atoms);
            if (g->match(ace_graph, perm)) {
                //printf("matched ACE\n");
                return true;
            }
        }
        if (h==nme_hash) {
            if (!g) g = Graph::create(mol, atoms);
            if (g->match(nme_graph, perm)) {
                //printf("matched NME\n");
                return true;
            }
        }
        return false;
    }

    NameMap types = {
          { "CA", AtomProBack }
        , { "C",  AtomProBack }
        , { "O",  AtomProBack }
        , { "N",  AtomProBack }
        , { "P",  AtomNucBack }
        , { "O1P",  AtomNucBack }
        , { "O2P",  AtomNucBack }
        , { "OP1",  AtomNucBack }
        , { "OP2",  AtomNucBack }
        , { "C3*",  AtomNucBack }
        , { "C3'",  AtomNucBack }
        , { "O3*",  AtomNucBack }
        , { "O3'",  AtomNucBack }
        , { "C4*",  AtomNucBack }
        , { "C4'",  AtomNucBack }
        , { "C5*",  AtomNucBack }
        , { "C5'",  AtomNucBack }
        , { "O5*",  AtomNucBack }
        , { "O5'",  AtomNucBack }
    };

    NameMap terms = {
          { "OT1", AtomProBack }
        , { "OT2", AtomProBack }
        , { "OX1", AtomProBack }
        , { "O1",  AtomProBack }
        , { "O2",  AtomProBack }
        , { "H5T", AtomNucBack }
        , { "H3T", AtomNucBack }
    };

    void analyze_residue(SystemPtr self, Id res) {
        std::vector<IdPair> perm;

        /* clear structure */
        IdList const& atoms = self->atomsForResidue(res);
        for (Id i=0; i<atoms.size(); i++) self->atomFAST(atoms[i]).type=AtomOther;
        self->residueFAST(res).type = ResidueOther;

        /* check for water */
        if (has_water_residue_name(self->residueFAST(res).name) || (
            Graph::hash(self,atoms)==waterhash &&
            Graph::create(self,atoms)->match(watergraph,perm))) {
            self->residueFAST(res).type = ResidueWater;
            return;
        }

        /* need at least four atoms to determine protein or nucleic */
        if (atoms.size()<4) return;

        if (is_protein_capping_group(self, atoms)) {
            self->residueFAST(res).type = ResidueProtein;
            return;
        }

        int npro=0, nnuc=0;
        std::set<std::string> names;

        Id ca_atm = BadId;
        Id c_atm = BadId;
        Id n_atm = BadId;
        for (Id i=0; i<atoms.size(); i++) {
            Id id = atoms[i];
            const atom_t& atm = self->atomFAST(id);
            const std::string& aname = atm.name;
            if (!names.insert(aname).second) continue;
            /* check for nucleic or protein backbone */
            NameMap::const_iterator iter=types.find(aname);
            AtomType atype=AtomOther;
            if (iter!=types.end()) {
                atype=iter->second;
                if (atype==AtomProBack) {
                    if (aname=="CA") ca_atm = id;
                    else if (aname=="C") c_atm = id;
                    else if (aname=="N") n_atm = id;
                }
            } else {
                /* try terminal names */
                iter=terms.find(aname);
                if (iter!=terms.end()) {
                    /* must be bonded to atom of the same type */
                    IdList const& bonds = self->bondsForAtom(id);
                    for (Id j=0; j<bonds.size(); j++) {
                        Id o = self->bondFAST(bonds[j]).other(id);
                        AtomType otype = self->atomFAST(o).type;
                        if (otype==iter->second) {
                            atype=otype;
                            break;
                        }
                    }
                }
            }
            self->atomFAST(id).type = atype;
            if (atype==AtomProBack) ++npro;
            if (atype==AtomNucBack) ++nnuc;
        }
        ResidueType rtype=ResidueOther;
        if (npro>=4 &&
            ca_atm!=BadId && c_atm!=BadId && n_atm != BadId &&
            !bad(self->findBond(ca_atm,n_atm)) &&
            !bad(self->findBond(ca_atm,c_atm)) &&
             bad(self->findBond(c_atm,n_atm))) {
            rtype=ResidueProtein;
        } else if (nnuc>=4) {
            rtype=ResidueNucleic;
        } else for (Id i=0; i<atoms.size(); i++) {
            self->atomFAST(atoms[i]).type = AtomOther;
        }
        if (rtype==ResidueProtein) {
            find_sidechain(self.get(), res, ca_atm);
        }
        self->residueFAST(res).type=rtype;
    }
}

void desres::msys::Analyze(SystemPtr self) {

    self->updateFragids();
    for (auto it=self->residueBegin(), e=self->residueEnd(); it!=e; ++it) {
        analyze_residue(self, *it);
    }
}


static const double def_bond = 1.0;
static const double def_angle = 109.5 * M_PI/180;
static const double def_dihedral = -120.0 * M_PI/180;

void desres::msys::GuessHydrogenPositions(SystemPtr mol, IdList const& hatoms) {

    /* Randomize orientation of water molecules, deterministically. */
    srand48(1999);

    /* partition by root atom */
    typedef std::map<Id,IdList> RootMap;
    RootMap map;
    for (Id h : hatoms) {
        if (mol->bondCountForAtom(h)!=1) continue;
        map[mol->bondedAtoms(h).at(0)].push_back(h);
    }

    /* process each root atom */
    for (RootMap::iterator it=map.begin(); it!=map.end(); ++it) {
        const Id r = it->first;
        const atom_t& root = mol->atom(r);
        IdList& hlist = it->second;
        std::sort(hlist.begin(), hlist.end());

        /* get other atoms bonded to root. */
        IdList c;
        for (Id b : mol->bondedAtoms(r)) {
            if (!std::binary_search(hlist.begin(), hlist.end(), b)) {
                c.push_back(b);
            }
        }

        if (hlist.size()==2) {
            atom_t& hyd1 = mol->atom(hlist[0]);
            atom_t& hyd2 = mol->atom(hlist[1]);
            if (c.empty()) {
                /* water */

                /* sample a point on the unit sphere */
                double u1 = drand48();
                double u2 = drand48();
                double z = 2*u1 - 1.0;
                double phi = 2*M_PI*u2;
                double R = sqrt(1-z*z);

                /* orient and place the first hydrogen */
                hyd1.x = root.x + def_bond * R * cos(phi);
                hyd1.y = root.y + def_bond * R * sin(phi);
                hyd1.z = root.z + def_bond * z;

                /* arbitrary reference position for the second hydrogen */
                Float A[3] = {hyd1.x, hyd1.y + def_bond, hyd1.z};
                /* position with random rotation about O-H1 */
                apply_dihedral_geometry(hyd2.pos(),
                        A,
                        mol->atom(hlist[0]).pos(),
                        mol->atom(r).pos(),
                        def_bond, def_angle, 2*M_PI*drand48());

            } else if (c.size()==1) {
                Id C = c[0];
                Id A = C;
                if (mol->bondCountForAtom(C)>1) {
                    for (auto id : mol->bondedAtoms(C)) {
                        A=id;
                        if (A!=r) break;
                    }
                }
                apply_dihedral_geometry(hyd1.pos(),
                        mol->atom(A).pos(),
                        mol->atom(C).pos(),
                        mol->atom(r).pos(),
                        def_bond, def_angle, -M_PI/2+def_dihedral);

                apply_dihedral_geometry(hyd2.pos(),
                        mol->atom(A).pos(),
                        mol->atom(C).pos(),
                        mol->atom(r).pos(),
                        def_bond, def_angle, -M_PI/2-def_dihedral);
            } else if (c.size()==2) {
                apply_dihedral_geometry(hyd1.pos(),
                        mol->atom(c[0]).pos(),
                        mol->atom(c[1]).pos(),
                        mol->atom(r).pos(),
                        def_bond, def_angle, def_dihedral);
                apply_dihedral_geometry(hyd2.pos(),
                        mol->atom(c[1]).pos(),
                        mol->atom(c[0]).pos(),
                        mol->atom(r).pos(),
                        def_bond, def_angle, def_dihedral);
            }

        } else if (hlist.size()==3) {
            if (c.size()>0) {
                Id C = c[0];
                Id A = C;
                if (mol->bondCountForAtom(C)>1) {
                    for (auto id : mol->bondedAtoms(C)) {
                        A=id;
                        if (A!=r) break;
                    }
                }
                apply_dihedral_geometry(mol->atom(hlist[0]).pos(),
                        mol->atom(A).pos(),
                        mol->atom(C).pos(),
                        mol->atom(r).pos(),
                        def_bond, def_angle, def_dihedral/2);
                apply_dihedral_geometry(mol->atom(hlist[1]).pos(),
                        mol->atom(hlist[0]).pos(),
                        mol->atom(C).pos(),
                        mol->atom(r).pos(),
                        def_bond, def_angle, def_dihedral);
                apply_dihedral_geometry(mol->atom(hlist[2]).pos(),
                        mol->atom(C).pos(),
                        mol->atom(hlist[0]).pos(),
                        mol->atom(r).pos(),
                        def_bond, def_angle, def_dihedral);

            }

        } else if (hlist.size()==1) {
            atom_t& hyd1 = mol->atom(hlist[0]);

            if (c.empty()) {
                hyd1.x = root.x;
                hyd1.y = root.y;
                hyd1.z = root.z + def_bond;

            } else if (c.size()==1) {
                /* find another atom bonded to c to define the plane */
                Id C = c[0];
                Id A = C;
                if (mol->bondCountForAtom(C)>1) {
                    for (auto id : mol->bondedAtoms(C)) {
                        A=id;
                        if (A!=r) break;
                    }
                }
                apply_dihedral_geometry(hyd1.pos(),
                        mol->atom(A).pos(),
                        mol->atom(C).pos(),
                        mol->atom(r).pos(),
                        def_bond, def_angle, 0);

            } else {
                Float ux=0, uy=0, uz=0;
                for (Id b : c) {
                    ux += mol->atom(r).x - mol->atom(b).x;
                    uy += mol->atom(r).y - mol->atom(b).y;
                    uz += mol->atom(r).z - mol->atom(b).z;
                }
                Float len = sqrt(ux*ux + uy*uy + uz*uz);
                if (len) {
                    ux /= len;
                    uy /= len;
                    uz /= len;
                }
                hyd1.x = root.x + def_bond*ux;
                hyd1.y = root.y + def_bond*uy;
                hyd1.z = root.z + def_bond*uz;
            }
        }
    }
}

// This is here rather than in term_table.cxx because term_table.cxx gets compiled into
// msys_core, and msys_core isn't versioned, so if you were to load an older msys core
// first you would hit a missing symbol error.
TermTablePtr desres::msys::ReplaceTableWithSortedTerms(TermTablePtr src) {

    // construct new table with temporary name
    auto sys = src->system();
    std::string name = src->name();
    std::string tmpname = name + ".tmp";
    while (sys->table(tmpname)) tmpname += "X";
    auto dst = sys->addTable(tmpname, src->atomCount(), src->params());
    // copy basic properties
    dst->category = src->category;
    dst->tableProps() = src->tableProps();
    // same atom ids
    IdList amap(sys->maxAtomId());
    std::iota(amap.begin(), amap.end(), 0);
    // same param ids
    IdList pmap(src->params()->paramCount());
    std::iota(pmap.begin(), pmap.end(), 0);
    // list of terms, sorted by atom ids
    IdList terms = src->terms();
    std::sort(terms.begin(), terms.end(),
            [&src](Id a, Id b) {
                auto aids = src->atoms(a);
                auto bids = src->atoms(b);
                if (aids == bids) {
                    Id pa = src->param(a);
                    Id pb = src->param(b);
                    return (src->params()->compare(pa, pb) < 0);
                }
                return (aids < bids);
            });


    AppendTerms(dst, src, amap, terms, pmap);
    // swap the tables
    src->destroy();
    dst->rename(name);

    return dst;
}
void desres::msys::ReplaceTablesWithSortedTerms(SystemPtr mol) {
    for (auto name : mol->tableNames()) {
        ReplaceTableWithSortedTerms(mol->table(name));
    }
}

