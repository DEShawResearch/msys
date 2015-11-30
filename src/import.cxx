#include "import.hxx"
#include <boost/algorithm/string.hpp> /* for boost::trim */

using namespace desres::msys;

void SystemImporter::initialize(IdList const& atoms) {
    resmap.clear();
    chnmap.clear();
    IdList reslist, chnlist;
    BOOST_FOREACH(Id atm, atoms) {
        Id res = sys->atom(atm).residue;
        Id chn = sys->residue(res).chain;
        reslist.push_back(res);
        chnlist.push_back(chn);
    }
    sort_unique(reslist);
    sort_unique(chnlist);

    BOOST_FOREACH(Id chn, chnlist) {
        chain_t const& chain = sys->chain(chn);
        chnmap[ChnKey(chain.ct, chain.name, chain.segid)] = chn;
        BOOST_FOREACH(Id res, sys->residuesForChain(chn)) {
            if (std::binary_search(reslist.begin(), reslist.end(), res)) {
                residue_t const& residue = sys->residue(res);
                resmap[ResKey(chn, residue.resid, residue.name, residue.insertion)] = res;
            }
        }
    }
}

bool SystemImporter::terminateChain(std::string chain, std::string segid, 
                                    Id ct ) {
    boost::trim(chain);
    boost::trim(segid);
    ChnMap::iterator it=chnmap.find(ChnKey(ct, chain, segid));
    if (it==chnmap.end()) return false;
    chnmap.erase(it);
    return true;
}

Id SystemImporter::addAtom(std::string chain, std::string segid,
                           int resnum, std::string resname,
                           std::string aname,
                           std::string insertion,
                           Id ct) {

    boost::trim(chain);
    boost::trim(segid);
    boost::trim(resname);
    boost::trim(insertion);
    boost::trim(aname);

    if (bad(ct)) MSYS_FAIL("Got ct=BadId");

    /* start a new ct if necessary */
    while (!sys->hasCt(ct)) {
        sys->addCt();
    }

    /* start a new chain if necessary */
    std::pair<ChnMap::iterator,bool> cp;
    cp = chnmap.insert(std::make_pair(ChnKey(ct,chain,segid),chnid));
    if (cp.second) {
        /* new chain/segid pair found, so start a new chain */
        chnid = sys->addChain(ct);
        sys->chainFAST(chnid).name = chain;
        sys->chainFAST(chnid).segid = segid;
        resid = BadId;
        cp.first->second = chnid;
    } else {
        /* use existing chain */
        chnid = cp.first->second;
    }

    /* start a new residue if necessary */
    std::pair<ResMap::iterator,bool> p;
    p = resmap.insert(std::make_pair(ResKey(chnid,resnum,resname,insertion), 
                resid));
    if (p.second) {
        /* new resname/resnum in this chain, so start a new residue. */
        resid = sys->addResidue(chnid);
        sys->residueFAST(resid).resid = resnum;
        sys->residueFAST(resid).name = resname;
        sys->residueFAST(resid).insertion = insertion;
        p.first->second = resid;
    } else {
        /* use existing residue */
        resid = p.first->second;
    }

    /* add the atom */
    Id atm = sys->addAtom(resid);
    sys->atomFAST(atm).name = aname;
    return atm;
}

