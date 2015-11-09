#ifndef desres_msys_atomsel_hxx
#define desres_msys_atomsel_hxx

#include "system.hxx"

namespace desres { namespace msys { 

    /* evaluate a vmd atom selection, returning the selected atoms. */
    IdList Atomselect(SystemPtr sys, const std::string& sel);

    /* evaluate with supplied positions and cell */
    IdList Atomselect(SystemPtr sys, const std::string& sel,
                      const float* pos, const double* cell);

}}

#endif
