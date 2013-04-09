#include "atomsel.hxx"

#ifndef MSYS_WITHOUT_ATOMSEL
#include "atomsel/vmd.hxx"
#include "atomsel/msys_keyword.hxx"
#endif

namespace desres { namespace msys { 

    IdList Atomselect(SystemPtr ptr, const std::string& txt) {
#ifndef MSYS_WITHOUT_ATOMSEL
        atomsel::VMD vmd(ptr);
        atomsel::PredicatePtr result = vmd.parse(txt);
        if (!result) {
            MSYS_FAIL("invalid selection:" << vmd.error);
        }
        atomsel::Selection s(atomsel::full_selection(ptr));
        result->eval(s);
        return s.ids();
#else
        MSYS_FAIL("No atom selection support");
#endif
    }

}}
