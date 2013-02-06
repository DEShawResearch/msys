#ifndef desres_msys_aromatic_new_hxx
#define desres_msys_aromatic_new_hxx

#include "system.hxx"

namespace desres { namespace msys {

    bool IsAromaticBond(SystemPtr sys, Id bond);
    bool IsAromaticAtom(SystemPtr sys, Id atom);

}}

#endif
