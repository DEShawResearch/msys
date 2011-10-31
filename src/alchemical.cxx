#include "alchemical.hxx"
#include "clone.hxx"
#include <sstream>
#include <stdexcept>
#include <boost/foreach.hpp>

using namespace desres::msys;

#define ERROR(x) do { \
    std::stringstream ss; \
    ss << "ERROR: " << x; \
    throw std::runtime_error(ss.str()); \
} while (0)


/* construct and return a mapping from atoms in the given input system to
 * the merged alchemical system, assuming 0-based ids.  */
static IdList map_atoms( IdList const& atoms, Id maxid ) {
    std::vector<Id> reals;
    IdList map;
    Id n=0;
    BOOST_FOREACH(Id id, atoms) {
        if (bad(id)) {
            ++n;
        } else if (id>=maxid) {
            ERROR("Invalid id " << id << " in atom map");
        } else {
            reals.push_back(id);
            map.push_back(map.size()+n);
        }
    }
    std::sort(reals.begin(), reals.end());
    reals.resize(std::unique(reals.begin(), reals.end())-reals.begin());
    if (reals.size() + n != atoms.size()) {
        ERROR("atom map contains duplicates, or some atoms not mapped");
    }
    if (reals.size() != maxid) {
        ERROR("not all atoms in system are mapped");
    }
    return map;
}
                    

SystemPtr desres::msys::CreateAlchemical( SystemPtr A, IdList const& aids,
                                          SystemPtr B, IdList const& bids) {

    if (aids.size() != bids.size()) {
        ERROR("atom maps have different sizes: " << aids.size() << ", "
                    << bids.size());
    }

    IdList mapA, mapB;
    try {
        mapA = map_atoms(aids, A->maxAtomId());
    }
    catch (std::exception& e) {
        ERROR("Invalid map for A sites:\n" << e.what());
    }
    try {
        mapB =map_atoms(bids, B->maxAtomId());
    }
    catch (std::exception& e) {
        ERROR("Invalid map for B sites:\n" << e.what());
    }

    /* create T based on A; add dummies.  */
    SystemPtr T = Clone(A, A->atoms());
    Id nd = aids.size() - A->atomCount(); /* number of dummies */
    for (Id i=0; i<nd; i++) {
        /* TODO: add to existing residue if possible based on the identity
         * of the corresponding B atom */
        Id chain = T->addChain();
        Id res = T->addResidue(chain);
        T->addAtom(res);
    }

    /* permute T to C such that the dummies are intermingled as specified
     * in the original atom map. */
    SystemPtr C;
    {
        IdList t(aids);
        Id d=A->maxAtomId();
        BOOST_FOREACH( Id& id, t ) {
            if (bad(id)) id=d++;
        }
        C = Clone(T, t);
    }

    return C;
}
