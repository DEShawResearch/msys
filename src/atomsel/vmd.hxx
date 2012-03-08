#ifndef desres_msys_atomsel_vmd_hxx
#define desres_msys_atomsel_vmd_hxx

#include "predicate.hxx"
#include "../system.hxx"

namespace desres { namespace msys { namespace atomsel { namespace vmd {

    typedef std::vector<const char *> StrList;
    PredicatePtr parse(const std::string& sel, SystemPtr sys, StrList& prev);

}}}}

#endif
