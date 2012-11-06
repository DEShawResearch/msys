#ifndef desres_msys_aromatic_hxx
#define desres_msys_aromatic_hxx

#include "system.hxx"

namespace desres { namespace msys {

    /* Balaban, Chem Rev, 2004, 104, 2777-2812 
       Implementation of *simple* aromaticity detection
       Given ring atom types {X[2e], Y[1e], Z[0e]}, a ring of size M 
       should satisfy:

       X+Y+Z=M
       X+0.5Y=2n+1  -> aromatic
       X+0.5Y=2n    -> antiaromatic

       The following additional restrictions are *NOT* considered here:
       1) sum of aromaticity constants should be -200< k[A] <200
       2) adjacent atoms of same type (for X and Z) destabilize the ring 
    */

    struct AromaticAtom {
        enum Type {
            X_TYPE=0,
            Y_TYPE,
            YEXT_TYPE,
            Z_TYPE,
            INVALID
        };

        /* FIXME: I have no idea what these inputs mean */
        static Type Classify(int nb, int a0, int b0, int b1, int be);
    };

    struct AromaticRing {
        enum Type {
            ANTIAROMATIC=0,
            NONAROMATIC,
            AROMATIC
        };

        /* Classify according to aromatic atom types */
        static Type Classify(int nX, int nY, int nYe, int nZ);
    };



    /* Find the number of atoms posessing the given AromaticAtomClassification
     * in the given ring.  Return true if all atoms were able to be classified,
     * false if not. */
    bool ClassifyAromaticAtoms(SystemPtr mol, IdList const& atoms,
                               int& nx, int& ny, int& nyext, int& nz);

    /* Get the aromatic ring type of the given atoms. */
    AromaticRing::Type ClassifyAromaticRing(SystemPtr mol, IdList const& atms);

    /* Compute ring planarity. */
    double ComputeRingPlanarity(SystemPtr mol, IdList const& atoms);

}}

#endif

