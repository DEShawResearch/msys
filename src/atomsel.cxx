#include "atomsel.hxx"
#include "atomsel/vmd.hxx"
#include "atomsel/msys_keyword.hxx"

namespace desres { namespace msys { 

    IdList Atomselect(SystemPtr ptr, const std::string& txt) {
        atomsel::VMD vmd(ptr);
        atomsel::PredicatePtr result = vmd.parse(txt);
        if (!result) {
            MSYS_FAIL("invalid selection:" << vmd.error);
        }
        atomsel::Selection s(atomsel::full_selection(ptr));
        result->eval(s);
        return s.ids();
    }

}}
