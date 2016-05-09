#include "atomsel.hxx"
#include "atomsel/vmd.hxx"
#include "atomsel/msys_keyword.hxx"

namespace desres { namespace msys { 

    IdList Atomselect(SystemPtr ptr, const std::string& txt) {
        return Atomselect(ptr, txt, nullptr, nullptr);
    }

    IdList Atomselect(SystemPtr ptr, const std::string& txt,
                      const float* pos, const double* cell) {

        atomsel::VMD vmd(ptr, pos, cell);
        try {
            atomsel::PredicatePtr result = vmd.parse(txt);
            if (!result) {
                MSYS_FAIL("invalid selection:" << vmd.error);
            }
            atomsel::Selection s(atomsel::full_selection(ptr));
            result->eval(s);
            return s.ids();
        }
        catch (std::exception& e) {
            MSYS_FAIL("Failed on selection: " << txt);
        }
    }

}}
