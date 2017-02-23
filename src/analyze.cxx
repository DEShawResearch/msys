#include "analyze.hxx"
#include "analyze/bond_orders.hxx"
#include "elements.hxx"
#include "graph.hxx"
#include "geom.hxx"
#include "contacts.hxx"
#include <stdio.h>
#include <queue>


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

    void AssignBondOrderAndFormalCharge(SystemPtr mol) {
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
                AssignBondOrderAndFormalCharge(mol, fragments[frags[0]]);
                continue;
            }

            /* We have multiple fragments with the same formula.  */
            for (Id frag : frags) {
                graphs[frag] = Graph::create(mol, fragments[frag]);
            }
            std::vector<IdPair> perm;
            while (!frags.empty()) {
                AssignBondOrderAndFormalCharge(mol, fragments[frags[0]]);
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
                                        int total_charge) {
#ifdef MSYS_WITHOUT_LPSOLVE
        MSYS_FAIL("LPSOLVE functionality was not included.");
#else
        if (atoms.empty()) return;
        BondOrderAssigner boa(mol, atoms);
        if (total_charge != INT_MAX) {
            boa.setTotalCharge(total_charge);
        }
        boa.solveIntegerLinearProgram();
        boa.assignSolutionToAtoms();
#endif
    }


    void GuessBondConnectivity(SystemPtr mol) {
        std::vector<Float> pos(3*mol->maxAtomId());
        static const double cutoff = 4.0;
        if (pos.empty()) return;
        IdList atoms(mol->atoms());
        for (Id i : atoms) {
            atom_t const& atom = mol->atom(i);
            pos[3*i  ] = atom.x;
            pos[3*i+1] = atom.y;
            pos[3*i+2] = atom.z;
        }
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

    IdList FindDistinctFragments(SystemPtr mol, MultiIdList const& fragments) {
        IdList result;
        /* will compute graphs lazily */
        std::vector<GraphPtr> graphs(fragments.size());
        typedef std::map<std::string, IdList> FragmentHash;
        FragmentHash fragment_hash;
        for (Id i=0; i<fragments.size(); i++) {
            fragment_hash[Graph::hash(mol, fragments[i])].push_back(i);
        }
        FragmentHash::iterator it;
        for (it=fragment_hash.begin(); it!=fragment_hash.end(); ++it) {
            /* unique formula -> unique fragment */
            IdList& frags = it->second;
            if (frags.size()==1) {
                result.push_back(frags[0]);
                continue;
            }
            /* must do isomorphism checks. */
            for (Id frag : frags) {
                graphs[frag] = Graph::create(mol, fragments[frag]);
            }
            std::vector<IdPair> perm;
            while (!frags.empty()) {
                result.push_back(frags[0]);
                IdList unmatched;
                GraphPtr ref = graphs[frags[0]];
                for (Id i=1; i<frags.size(); i++) {
                    GraphPtr sel = graphs[frags[i]];
                    if (!ref->match(sel, perm)) {
                        unmatched.push_back(frags[i]);
                    }
                }
                frags.swap(unmatched);
            }
        }
        std::sort(result.begin(), result.end());
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
    struct _ { 
        _() {
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
    } static_initializer;

    typedef std::map<std::string,AtomType> NameMap;

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

