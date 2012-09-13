#include "alchemical.hxx"
#include "clone.hxx"
#include "schema.hxx"

#include <stdio.h>
#include <stack>
#include <math.h>

using namespace desres::msys;

namespace {

    Id copy_param(ParamTablePtr dstparams, 
                  ParamTablePtr srcparams,
                  Id srcid) {

        if (dstparams == srcparams) {
            return srcparams->duplicate(srcid);
        }
        Id dstid = dstparams->addParam();
        for (Id dstindex=0; dstindex<dstparams->propCount(); dstindex++) {
            String prop = dstparams->propName(dstindex);
            Id srcindex = srcparams->propIndex(prop);
            if (bad(srcindex)) continue;
            dstparams->value(dstid, dstindex)=srcparams->value(srcid,srcindex);
        }
        return dstid;
    }

    typedef std::pair<int,int> IntPair;
    typedef std::map<IdList, IntPair> BlockMap;

    TermTablePtr copy_table_template(SystemPtr mol, TermTablePtr B) {
        TermTablePtr A = mol->addTable(B->name(), B->atomCount());
        ParamTablePtr ap = A->params();
        ParamTablePtr bp = B->params();
        for (Id i=0; i<bp->propCount(); i++) {
            ap->addProp(bp->propName(i), bp->propType(i));
        }
        for (Id i=0; i<B->termPropCount(); i++) {
            A->addTermProp(B->termPropName(i), B->termPropType(i));
        }
        return A;
    }

    void canonicalize(IdList& ids) {
        size_t n = ids.size();
        Id a=0,b=0;
        switch (n) {
            case 2: 
                a=ids[0];
                b=ids[1];
                break;
            case 3:
                a=ids[0];
                b=ids[2];
                break;
            case 4:
                a=ids[1];
                b=ids[2];
                break;
            default:
                ;
        }
        if (a>b) std::reverse(ids.begin(), ids.end());
    }

    BlockMap make_pairs(SystemPtr mol, 
                         TermTablePtr atable, TermTablePtr btable,
                         IdList const& alist,
                         IdList const& blist,
                         IdList const& b2a,
                         IdList const& a2b) {

        BlockMap M;
        if (!btable) return M;
        if (!atable) atable = copy_table_template(mol, btable);

        /* index the B terms */
        IdList bterms = btable->findWithAny(blist);
        BOOST_FOREACH(Id i, bterms) {
            IdList ids = btable->atoms(i);
            BOOST_FOREACH(Id& id, ids) id = b2a[id];
            canonicalize(ids);
            M[ids]=IntPair(0,i+1);
        }

        /* index the A terms, combining them with any existing B terms */
        IdList aterms = atable->findWithAny(alist);
        BOOST_FOREACH(Id i, aterms) {
            IdList ids = atable->atoms(i);
            canonicalize(ids);
            BlockMap::iterator p = M.insert(
                    std::make_pair(ids,IntPair(i+1,0))).first;
            p->second.first = i+1;
        }

        return M;
    }

    Id make_full_pair(TermTablePtr pairs, Id ai, Id aj) {
        SystemPtr mol = pairs->system();
        TermTablePtr nb = mol->table("nonbonded");
        Id ti = nb->findWithAll(IdList(1,ai)).at(0);
        Id tj = nb->findWithAll(IdList(1,aj)).at(0);
        Float qi = mol->atom(ai).charge;
        Float qj = mol->atom(aj).charge;
        Id param = pairs->params()->addParam();
        if (pairs->name()=="pair_12_6_es") {
            if (mol->nonbonded_info.vdw_funct != "vdw_12_6") {
                MSYS_FAIL("Expected vdw_funct " << mol->nonbonded_info.vdw_funct << " with pairs_12_6_es table");
            }
            Float si = nb->propValue(ti,"sigma");
            Float sj = nb->propValue(tj,"sigma");
            Float ei = nb->propValue(ti,"epsilon");
            Float ej = nb->propValue(tj,"epsilon");
            Float sij=0, eij=0;
            if (mol->nonbonded_info.vdw_rule == "geometric") {
                sij = sqrt(si*sj);
                eij = sqrt(ei*ej);
            } else if (mol->nonbonded_info.vdw_rule == "arithmetic/geometric") {
                sij = 0.5*(si+sj);
                eij = sqrt(ei*ej);
            } else {
                MSYS_FAIL("Unsupported vdw_rule " << mol->nonbonded_info.vdw_rule);
            }
            Float aij = pow(sij,12) * eij * 4.0;
            Float bij = pow(sij, 6) * eij * 4.0;
            Float qij = qi * qj;
            pairs->params()->value(param,"aij") = aij;
            pairs->params()->value(param,"bij") = bij;
            pairs->params()->value(param,"qij") = qij;
        } else {
            MSYS_FAIL("Unsupported pairs table type " << pairs->name());
        }
        return param;
    }

    BlockMap make_exclmap(SystemPtr mol, 
                         TermTablePtr atable, TermTablePtr btable,
                         IdList const& alist,
                         IdList const& blist,
                         IdList const& b2a,
                         IdList const& a2b) {

        BlockMap M;
        if (!btable) return M;
        std::string ptype = "pair_12_6_es";
        TermTablePtr pairs = AddTable(mol,ptype);
        TermTablePtr pairsB = AddTable(btable->system(), ptype);
        TermTablePtr alcpairs = AddTable(mol, "alchemical_"+ptype);

        if (!atable) atable = AddTable(mol,"exclusion");

        /* index the B terms */
        IdList bterms = btable ? btable->findWithAny(blist) : IdList();
        BOOST_FOREACH(Id i, bterms) {
            IdList ids = btable->atoms(i);
            BOOST_FOREACH(Id& id, ids) id = b2a[id];
            canonicalize(ids);
            M[ids]=IntPair(-1,i+1);
        }

        /* index the A terms combining with existing B */
        IdList aterms = atable->findWithAny(alist);
        BOOST_FOREACH(Id i, aterms) {
            IdList ids = atable->atoms(i);
            canonicalize(ids);
            BlockMap::iterator p = M.insert(
                    std::make_pair(ids,IntPair(i+1,-1))).first;
            p->second.first = i+1;
        }

        /* process merged terms */
        for (BlockMap::iterator it=M.begin(); it!=M.end(); ++it) {
            IdList const& ids = it->first;
            Id ai = ids[0];
            Id aj = ids[1];
            int& ti = it->second.first;
            int& tj = it->second.second;

            if (tj==-1) {
                Id bi = a2b[ai];
                Id bj = a2b[aj];
                /* An excluded pair of atoms in A mapped onto an unexcluded
                 * pair of atoms in B.  */
                if (std::binary_search(alist.begin(), alist.end(), ai) &&
                    std::binary_search(alist.begin(), alist.end(), aj) &&
                    !bad(bi) && !bad(bj)) {
                    //printf("offset exclusion %u-%u in A by pair in B: %u-%u\n",
                            //ai,aj, bi,bj);
                    /* remove non-alchemical pair if it exists */
                    BOOST_FOREACH(Id id, pairs->findWithAll(ids)) {
                        pairs->delTerm(id);
                    }
                    /* find or create the alchemical pair */
                    Id param = BadId;
                    BOOST_FOREACH(Id id, alcpairs->findWithAll(ids)) {
                        param = alcpairs->param(id);
                        break;
                    }
                    if (bad(param)) {
                        param = alcpairs->params()->addParam();
                        alcpairs->addTerm(ids, param);
                    }
                    /* assign alchemical state */
                    Id full_pair = make_full_pair(pairsB, bi, bj);
                    for (Id i=0; i<pairsB->params()->propCount(); i++) {
                        std::string prop = pairsB->params()->propName(i);
                        prop += "B";
                        alcpairs->params()->value(param,prop) = 
                           pairsB->params()->value(full_pair,i);
                    }
                }
            } else if (ti==-1) {
                Id bi = a2b[ai];
                Id bj = a2b[aj];
                //printf("B: consider B %u %u mapped to A %u %u\n",
                        //bi,bj,ai,aj);
                /* An excluded pair of atoms in B mapped onto an unexcluded
                 * pair of atoms in A.  */
                /* FIXME: lots of code duplication with the first path */
                if (std::binary_search(alist.begin(), alist.end(), ai) &&
                    std::binary_search(alist.begin(), alist.end(), aj) &&
                    !bad(bi) && !bad(bj)) {
                    //printf("offset exclusion %u-%u in B by pair in A\n", ai,aj);
                    /* remove non-alchemical pair if it exists */
                    BOOST_FOREACH(Id id, pairs->findWithAll(ids)) {
                        pairs->delTerm(id);
                    }
                    /* find or create the alchemical pair */
                    Id param = BadId;
                    BOOST_FOREACH(Id id, alcpairs->findWithAll(ids)) {
                        param = alcpairs->param(id);
                        break;
                    }
                    if (bad(param)) {
                        param = alcpairs->params()->addParam();
                        alcpairs->addTerm(ids, param);
                    }
                    /* assign alchemical state */
                    Id full_pair = make_full_pair(pairs, ai, aj);
                    for (Id i=0; i<pairs->params()->propCount(); i++) {
                        std::string prop = pairs->params()->propName(i);
                        prop += "A";
                        alcpairs->params()->value(param,prop) = 
                           pairs->params()->value(full_pair,i);
                    }
                }
            }
        }

        /* add more exclusions based on distance cutoff */
        const Float cut = 6.0;
        const Float r2 = cut*cut;
        SystemPtr bmol = btable->system();
        BOOST_FOREACH(Id ai, alist) {
            if (!bad(a2b.at(ai))) continue;
            atom_t const& a = mol->atom(ai);
            BOOST_FOREACH(Id bi, blist) {
                if (std::binary_search(alist.begin(),alist.end(),b2a.at(bi))) {
                    continue;
                }
                atom_t const& b = bmol->atom(bi);
                Float dx = b.x-a.x;
                Float dy = b.y-a.y;
                Float dz = b.z-a.z;
                Float d2 = dx*dx + dy*dy + dz*dz;
                if (d2<r2) {
                    IdList ids(2);
                    ids[0] = ai;
                    ids[1] = b2a.at(bi);
                    M[ids] = IntPair(-1,-1);
                }
            }
        }
        return M;
    }

    BlockMap make_block(SystemPtr mol, 
                         TermTablePtr atable, TermTablePtr btable,
                         IdList const& alist,
                         IdList const& blist,
                         IdList const& b2a,
                         IdList const& a2b) {

        BlockMap M;
        if (!btable) return M;
        if (!atable) atable = copy_table_template(mol, btable);

        /* index the B terms */
        IdList bterms = btable->findWithAny(blist);
        BOOST_FOREACH(Id i, bterms) {
            IdList ids = btable->atoms(i);
            BOOST_FOREACH(Id& id, ids) id = b2a[id];
            canonicalize(ids);
            M[ids]=IntPair(-1,i+1);
        }

        /* index the A terms, combining them with any existing B terms */
        IdList aterms = atable->findWithAny(alist);
        BOOST_FOREACH(Id i, aterms) {
            IdList ids = atable->atoms(i);
            canonicalize(ids);
            BlockMap::iterator p = M.insert(
                    std::make_pair(ids,IntPair(i+1,-1))).first;
            p->second.first = i+1;
        }

        /* Change ti,tj=0 to -1 for "kept" terms */
        for (BlockMap::iterator it=M.begin(); it!=M.end(); ++it) {
            IdList const& ids = it->first;
            int& ti = it->second.first;
            int& tj = it->second.second;

            if (ti==-1) {
                BOOST_FOREACH(Id const& i, ids) { 
                    if (std::binary_search(alist.begin(),alist.end(),i)) {
                        ti=0;
                        break;
                    }
                }
            } else if (tj==-1) {
                BOOST_FOREACH(Id const& id, ids) {
                    if (!bad(a2b.at(id))) {
                        tj=0;
                        break;
                    }
                }
            }
        }
        return M;
    }

    Id copy_alchemical(ParamTablePtr dst, 
                       ParamTablePtr srcA, Id idA,
                       ParamTablePtr srcB, Id idB) {
        Id id = dst->addParam();
        if (!bad(idA)) {
            for (Id i=0; i<srcA->propCount(); i++) {
                std::string prop = srcA->propName(i);
                Id col = dst->addProp(prop+'A', srcA->propType(i));
                dst->value(id,col) = srcA->value(idA,i);
            }
        }
        if (!bad(idB)) {
            for (Id i=0; i<srcB->propCount(); i++) {
                std::string prop = srcB->propName(i);
                Id col = dst->addProp(prop+'B', srcB->propType(i));
                dst->value(id,col) = srcB->value(idB,i);
            }
        }
        return id;
    }

    /* true if, for a property x in params, there are properties xA, xB
     * in alc such that value(p, xA)!=value(p, xB) */
    bool alchemically_different(ParamTablePtr alc, Id p, ParamTablePtr src) {
        for (Id i=0; i<src->propCount(); i++) {
            std::string prop = src->propName(i);
            Id colA = alc->propIndex(prop+"A");
            Id colB = alc->propIndex(prop+"B");
            if (!bad(colA) && !bad(colB)) {
                if (alc->value(p,colA) != alc->value(p,colB)) {
                    return true;
                }
            }
        }
        return false;
    }

    void make_alchemical_excl(TermTablePtr table, BlockMap const& block) {
        for (BlockMap::const_iterator it=block.begin(); it!=block.end(); ++it) {
            IdList const& ids = it->first;
            int const& ti = it->second.first;
            int const& tj = it->second.second;
            if (ti==-1 && (tj>0 || tj==-1)) {
                table->addTerm(ids, BadId);
            }
        }
    }

    void make_alchemical(TermTablePtr atable, TermTablePtr btable,
                         BlockMap const& block, 
                         bool avoid_noops, std::string const& keeper="") {

        if (!btable) return;
        Id b_constrained = btable->termPropIndex("constrained");
        Id a_constrained = BadId;
        if (!bad(b_constrained)) {
            a_constrained = atable->addTermProp("constrained", IntType);
        }
        std::string alcname = "alchemical_" + btable->name();
        ParamTablePtr params = atable->params();
        TermTablePtr alc = AddTable(atable->system(), alcname);
        ParamTablePtr paramsB = alc->params();

        for (BlockMap::const_iterator it=block.begin(); it!=block.end(); ++it) {
            IdList const& ids = it->first;
            int ta = it->second.first;
            int tb = it->second.second;

            if (ta>0 && tb==-1) {
                /* Keep the A state intact, no alchemical transformation */

            } else if (ta>0 && tb==0) {
                /* disappear the B state */
                Id p = copy_alchemical(paramsB, 
                                       params, atable->param(ta-1),
                                       ParamTablePtr(), BadId);
                if (!keeper.empty()) {
                    Id colA = paramsB->propIndex(keeper+'A');
                    Id colB = paramsB->propIndex(keeper+'B');
                    paramsB->value(p,colB)=paramsB->value(p,colA);
                }
                alc->addTerm(atable->atoms(ta-1), p);
                atable->delTerm(ta-1);

            } else if (ta==-1 && tb>0) {
                /* copy parameters from state B */
                Id p = copy_param(params, btable->params(), btable->param(tb-1));
                Id term = atable->addTerm(ids, p);
                if (!bad(a_constrained)) {
                    atable->termPropValue(term,a_constrained) =
                        btable->termPropValue(tb-1, b_constrained);
                }

            } else if (ta==0 && tb>0) {
                /* disappear the A state */
                Id p = copy_alchemical(paramsB, 
                                       ParamTablePtr(), BadId,
                                       btable->params(), btable->param(tb-1));
                if (!keeper.empty()) {
                    Id colA = paramsB->propIndex(keeper+'A');
                    Id colB = paramsB->propIndex(keeper+'B');
                    paramsB->value(p,colA)=paramsB->value(p,colB);
                }
                alc->addTerm(ids, p);

            } else if (ta>0 && tb>0) {
                /* A state morphs into B state */
                Id p = copy_alchemical(paramsB,
                                       atable->params(), atable->param(ta-1),
                                       btable->params(), btable->param(tb-1));
                if (!avoid_noops || 
                        alchemically_different(paramsB, p, atable->params())) {
                    alc->addTerm(ids, p);
                    atable->delTerm(ta-1);
                }

            } else {
                MSYS_FAIL("Unsupported mapping in " << atable->name()
                        << ": ta=" << ta << ", tb=" << tb);
            }
        }
    }

    IdList two_neighbors(Id i, SystemPtr mol, IdSet const& mapped,
            IdList const& idmap) {
        /* Return a list of 1, 2 or 3 atoms [i,j,k] such that i is 
         * bonded to j, j is bonded to k, and k!=i.  */
        IdList ret(1,idmap.at(i));
        bool type = mapped.count(i);

        /* find neighbors of j of the same type, along with their degree */
        IdList nbrs, degs, bonded = mol->bondedAtoms(i);
        std::sort(bonded.begin(), bonded.end());
        BOOST_FOREACH(Id j, bonded) {
            if (mapped.count(j)!=type) continue;
            Id deg=0;
            BOOST_FOREACH(Id k, mol->bondedAtoms(j)) {
                if (mapped.count(k)!=type) continue;
                ++deg;
            }
            nbrs.push_back(j);
            degs.push_back(deg);
        }
        if (nbrs.empty()) return ret;

        //printf("    neighbors of %u: ", idmap.at(i));
        //for (Id j=0; j<nbrs.size(); j++) printf("[%u, %u] ", nbrs[j], degs[j]);
        //printf("\n");

        /* keep the one with the highest degree */
        Id ind = std::max_element(degs.begin(), degs.end())-degs.begin();
        Id j = idmap.at(nbrs[ind]);
        ret.push_back(j);
        //printf("    got position %u with max degree, atom %u\n", ind,j);

        /* now try to find the second atom */
        BOOST_FOREACH(Id k, mol->bondedAtoms(nbrs[ind])) {
            if (mapped.count(i)==type && k!=i) {
                ret.push_back(idmap.at(k));
                //printf("        adding second atom %u\n", idmap.at(k));
                break;
            }
        }
        return ret;
    }

    void keep(BlockMap& map, Id ai, Id aj=BadId, Id ak=BadId, Id al=BadId) {
        IdList key(1,ai);
        if (!bad(aj)) key.push_back(aj);
        if (!bad(ak)) key.push_back(ak);
        if (!bad(al)) key.push_back(al);
        canonicalize(key);
        BlockMap::iterator it=map.find(key);
        if (it==map.end()) return;
        int& ti = it->second.first;
        int& tj = it->second.second;
        if (ti==0) ti=-1;
        else if (tj==0) tj=-1;
        else MSYS_FAIL("Invalid term: [" << ti << ", " << tj << "]");
    }

    void find_connected_dummies(Id i, SystemPtr mol, IdSet const& mapped,
                                IdSet& connected) {
        /* Find all the dummy atoms that are connected to i - directly or 
         * indirectly, through covalent bonds between dummy atoms, and set 
         * the corresponding connected bits to true for these dummy atoms.
         */
        std::stack<Id> S;
        S.push(i);
        std::vector<bool> visited(mol->maxAtomId(), false);
        while (!S.empty()) {
            Id k = S.top();
            S.pop();
            connected.insert(k);
            visited.at(k)=true;
            BOOST_FOREACH(Id j, mol->bondedAtoms(k)) {
                if (visited.at(j)) continue;
                if (!mapped.count(j)) continue;
                S.push(j);
            }
        }

    }

    void find_kept_atoms(SystemPtr A, IdList const& alist, IdSet const& mappedA,
                         IdList const& idmap,
                         BlockMap& bondmap,
                         BlockMap& anglmap,
                         BlockMap& dihemap) {

        IdSet connected;
        BOOST_FOREACH(Id i, alist) {
            if (mappedA.count(i) || connected.count(i)) continue;
            IdList jlist = A->bondedAtoms(i);
            /* FIXME: for comparison purposes only.  */
            std::sort(jlist.begin(), jlist.end());
            BOOST_FOREACH(Id j, jlist) {
                if (!mappedA.count(j)) continue;
                IdList n1 = two_neighbors(i,A,mappedA,idmap);
                IdList n2 = two_neighbors(j,A,mappedA,idmap);
                //printf("  using template [");
                //BOOST_FOREACH(Id id,n1) printf("%u ", id);
                //printf("] , [");
                //BOOST_FOREACH(Id id,n2) printf("%u ", id);
                //printf("]\n");

                if (1)            keep(bondmap, n1[0], n2[0]);
                if (n1.size()>=2) keep(anglmap, n1[1],n1[0],n2[0]);
                //if (n1.size()>=3) keep(dihemap, n1[2],n1[1],n1[0],n2[0]);
                if (n2.size()>=2) keep(anglmap, n2[1],n2[0],n1[0]);
                //if (n2.size()>=3) keep(dihemap, n2[2],n2[1],n2[0],n1[0]);
                break;
            }
            find_connected_dummies(i,A,mappedA,connected);
        }
    }

    void stage2(SystemPtr A, SystemPtr B,
                IdList const& alist,
                IdList const& blist,
                IdList const& b2a,
                IdList const& a2b,
                BlockMap& bondmap,
                BlockMap& anglmap,
                BlockMap& dihemap) {

        /* find which A atoms have real analogs in B, and vice versa */
        IdSet mappedA, mappedB;
        for (Id i=0; i<alist.size(); i++) {
            Id a = alist[i];
            Id b = a2b.at(a);
            if (!bad(b)) {
                mappedA.insert(a);
                mappedB.insert(b);
            }
        }
        IdList a2a(A->atomCount(), BadId);
        for (Id i=0; i<a2a.size(); i++) a2a[i]=i;
        //printf("stage A\n");
        find_kept_atoms(A,alist,mappedA, a2a, bondmap, anglmap, dihemap);
        //printf("stage B\n");
        find_kept_atoms(B,blist,mappedB, b2a, bondmap, anglmap, dihemap);
    }
}

namespace {
    void merge_constraints(SystemPtr A, SystemPtr B,
                           IdList const& alist, IdList const& blist,
                           IdList const& b2a, IdList const& a2b) {

        /* hash constraints by first atom */
        typedef std::pair<TermTablePtr, Id> term_t;
        typedef std::map<Id, term_t> TermMap;
        TermMap map;
        BOOST_FOREACH(String name, A->tableNames()) {
            TermTablePtr table = A->table(name);
            if (table->category!=CONSTRAINT) continue;
            BOOST_FOREACH(Id t, table->findWithAny(alist)) {
                Id a = table->atom(t,0);
                if (!map.insert(std::make_pair(a, term_t(table,t))).second) {
                    MSYS_FAIL("Overlapping constraint for atom " << a);
                }
            }
        }

        /* find constraints in B that overlap with A */
        BOOST_FOREACH(String name, B->tableNames()) {
            TermTablePtr btable = B->table(name);
            if (btable->category!=CONSTRAINT) continue;
            BOOST_FOREACH(Id t, btable->findWithAny(blist)) {
                Id a = b2a.at(btable->atom(t,0));
                TermMap::const_iterator it = map.find(a);
                if (it==map.end()) {
                    /* copy constraint from B to A */
                    TermTablePtr atable = AddTable(A, btable->name());
                    IdList ids(btable->atoms(t));
                    //printf("copy from %s: ids size %lu\n", 
                            //atable->name().c_str(), ids.size());
                    Id p = copy_param(atable->params(), 
                                      btable->params(), 
                                      btable->param(t));
                    BOOST_FOREACH(Id& id, ids) id=b2a.at(id);
                    atable->addTerm(ids, p);
                } else {
                    /* merge constraint from B into A */
                    TermTablePtr atable = it->second.first;
                    Id aterm = it->second.second;
                    IdList aids = atable->atoms(aterm);
                    IdList bids = btable->atoms(t);
                    IdList ids(aids), param_indices;
                    for (Id i=1; i<bids.size(); i++) {
                        Id id = bids[i];
                        id = b2a.at(id);
                        if (std::find(aids.begin(),aids.end(),id)==aids.end()) {
                            ids.push_back(id);
                            param_indices.push_back(i-1);
                        }
                    }
                    //printf("merge B %s: ", btable->name().c_str());
                    //BOOST_FOREACH(Id id, bids) printf("%u ", id);
                    //printf("into A %s: ", atable->name().c_str());
                    //BOOST_FOREACH(Id id, aids) printf("%u ", id);
                    //printf("\n  merged: ");
                    //BOOST_FOREACH(Id id, ids) printf("%u ", id);
                    //printf("\n");

                    std::stringstream ss;
                    ss << "constraint_ah" << ids.size()-1;
                    TermTablePtr c = AddTable(A,ss.str());
                    Id cparam = c->params()->addParam();
                    c->addTerm(ids, cparam);
                    ParamTablePtr bparams = btable->params();
                    for (Id i=0; i<atable->params()->propCount(); i++) {
                        c->params()->value(cparam,i) = atable->propValue(aterm,i);
                    }
                    for (Id i=0; i<param_indices.size(); i++) {
                        c->params()->value(cparam,i+atable->params()->propCount())
                            = btable->propValue(t,param_indices[i]);
                    }
                    atable->delTerm(aterm);
                }
            }
        }
    }
}

SystemPtr desres::msys::MakeAlchemical( SystemPtr A, SystemPtr B,
                                        std::vector<IdPair> pairs,
                                        bool avoid_noops) {

    /* We will be modifying A, so make a copy */
    A = Clone(A, A->atoms());

    /* Add custom atom properties from B */
    for (Id i=0; i<B->atomPropCount(); i++) {
        A->addAtomProp(B->atomPropName(i), B->atomPropType(i));
    }

    /* configure alchemical nonbonded table */
    TermTablePtr nbA = A->table("nonbonded");
    TermTablePtr nbB = B->table("nonbonded");
    TermTablePtr alc = A->addTable("alchemical_nonbonded", 1, nbA->params());
    alc->category = NONBONDED;
    alc->addTermProp("chargeB", FloatType);
    alc->addTermProp("moiety",  IntType);
    /* zero nonbonded for dummy atoms */
    Id zero = nbA->params()->addParam();

    /* construct mappings from B to A and from A to B, and the list of 
     * alchemical A and B atoms.  Add dummy atoms to A when B maps to BadId. */
    IdList a2b(A->atomCount(), BadId), b2a(B->atomCount(), BadId);
    IdList alist, blist, dummies;
    BOOST_FOREACH(IdPair const& p, pairs) {
        Id ai = p.first;
        Id bi = p.second;
        if (bad(ai) && bad(bi)) MSYS_FAIL("Cannot map dummy to dummy");
        if (!bad(ai)) alist.push_back(ai);
        if (!bad(bi)) blist.push_back(bi);
    }
    if (sort_unique(alist)) MSYS_FAIL("Duplicate A state ids in atom map");
    if (sort_unique(blist)) MSYS_FAIL("Duplicate B state ids in atom map");

    /* ensure that A and B contain no dangling terms; i.e. terms with 
     * both mapped and unmapped atoms. */
    BOOST_FOREACH(String name, A->tableNames()) {
        TermTablePtr table = A->table(name);
        BOOST_FOREACH(Id term, table->findWithAny(alist)) {
            BOOST_FOREACH(Id atom, table->atoms(term)) {
                if (!std::binary_search(alist.begin(), alist.end(), atom)) {
                    MSYS_FAIL("State A table " << name << " term " << term 
                            << " contains unmapped atom " << atom);
                }
            }
        }
    }
    BOOST_FOREACH(String name, B->tableNames()) {
        TermTablePtr table = B->table(name);
        BOOST_FOREACH(Id term, table->findWithAny(blist)) {
            BOOST_FOREACH(Id atom, table->atoms(term)) {
                if (!std::binary_search(blist.begin(), blist.end(), atom)) {
                    MSYS_FAIL("State B table " << name << " term " << term 
                            << " contains unmapped atom " << atom);
                }
            }
        }
    }

    SystemImporter imp(A);
    imp.initialize(alist);
    BOOST_FOREACH(IdPair const& p, pairs) {
        Id ai = p.first;
        Id bi = p.second;

        if (bad(ai)) {
            /* dummy atom for the A state */
            a2b.push_back(bi);
            atom_t const& batm = B->atom(bi);
            residue_t const& bres = B->residue(batm.residue);
            chain_t const& bchn = B->chain(bres.chain);
            ai = imp.addAtom(bchn.name, bchn.segid,
                             bres.resid, bres.name,
                             batm.name);
            dummies.push_back(ai);
            atom_t& atom = A->atom(ai);
            atom.mass = batm.mass;
            atom.x = batm.x;
            atom.y = batm.y;
            atom.z = batm.z;
            atom.vx = batm.vx;
            atom.vy = batm.vy;
            atom.vz = batm.vz;
            atom.atomic_number = batm.atomic_number;
            for (Id p=0; p<B->atomPropCount(); p++) {
                Id q = A->atomPropIndex(B->atomPropName(p));
                A->atomPropValue(ai,q) = B->atomPropValue(bi,p);
            }
            b2a[bi] = ai;
            nbA->addTerm(IdList(1,ai), zero);
            /* copy the nonbonded type from B for the alchemical state */
            Id p = copy_param(nbA->params(), nbB->params(), nbB->param(bi));
            Id alcid = alc->addTerm(IdList(1,ai), p);
            /* copy charge in B to chargeB */
            alc->termPropValue(alcid, 0)=batm.charge;
        } else {
            /* alchemical B state */
            if (bad(bi)) {
                alc->addTerm(IdList(1,ai), zero);
            } else {
                b2a[bi]=ai;
                a2b[ai]=bi;
                Id p = copy_param(nbA->params(), nbB->params(), nbB->param(bi));
                bool add_alc_term = true;
                if (avoid_noops) {
                    Id tA = nbA->findWithAll(IdList(1,ai)).at(0);
                    Id pA = nbA->param(tA);
                    if (A->atom(ai).charge == B->atom(bi).charge && 
                            !nbA->params()->compare(p,pA)) {
                        add_alc_term = false;
                    }
                }
                if (add_alc_term) {
                    Id alcid = alc->addTerm(IdList(1,ai), p);
                    alc->termPropValue(alcid, 0)=B->atom(bi).charge;
                }
            }
        }
    }

    /* handle bonded terms */
    BlockMap bondmap, anglmap, dihemap, imprmap, pairmap, exclmap;
    String ff;
    ff = "stretch_harm";
    bondmap=make_block(A,A->table(ff),B->table(ff),alist,blist,b2a,a2b);
    ff = "angle_harm";
    anglmap=make_block(A,A->table(ff),B->table(ff),alist,blist,b2a,a2b);
    ff = "dihedral_trig";
    dihemap=make_block(A,A->table(ff),B->table(ff),alist,blist,b2a,a2b);
    ff="improper_harm";
    imprmap=make_block(A,A->table(ff),B->table(ff),alist,blist,b2a,a2b);
    ff="pair_12_6_es";
    pairmap=make_pairs(A,A->table(ff),B->table(ff),alist,blist,b2a,a2b);

    stage2(A,B, alist, blist, b2a, a2b, bondmap, anglmap, dihemap);

    ff = "stretch_harm";
    make_alchemical(A->table(ff), B->table(ff), bondmap, avoid_noops, "r0");
    ff = "angle_harm";
    make_alchemical(A->table(ff), B->table(ff), anglmap, avoid_noops);
    ff = "dihedral_trig";
    make_alchemical(A->table(ff), B->table(ff), dihemap, avoid_noops);
    ff="improper_harm";
    make_alchemical(A->table(ff), B->table(ff), imprmap, avoid_noops);
    ff="pair_12_6_es";
    make_alchemical(A->table(ff), B->table(ff), pairmap, avoid_noops);

    /* handle bonds */
    BOOST_FOREACH(Id bi, blist) {
        Id ai = b2a.at(bi);
        assert(!bad(ai));
        BOOST_FOREACH(Id bnd, B->bondsForAtom(bi)) {
            bond_t const& bond = B->bond(bnd);
            Id bj = bond.other(bi);
            Id aj = b2a.at(bj);
            if (!bad(aj)) {
                A->addBond(ai,aj);
            }
        }
    }

    /* exclusions and pairs */
    ff = "exclusion";
    exclmap = make_exclmap(A,A->table(ff),B->table(ff),alist,blist,b2a,a2b);
    make_alchemical_excl(A->table(ff), exclmap);

    /* constraints */
    merge_constraints(A,B,alist,blist,b2a,a2b);

    /* coalesce, sort, and return */
    A->coalesceTables();

    /* copy the non-dummy atoms */
    IdList ids(A->atomCount()), atoms = A->atoms();
    ids.resize(
            std::set_difference(atoms.begin(), atoms.end(), 
                dummies.begin(), dummies.end(), ids.begin()) - ids.begin());
    /* insert the dummies after the last alist atom */
    IdList::iterator pos = std::find(ids.begin(), ids.end(), alist.back());
    ids.insert(pos+1, dummies.begin(), dummies.end());
    assert(ids.size()==atoms.size());
    return Clone(A, ids);
}

