#include "atomsel.hxx"

#ifndef MSYS_WITHOUT_ATOMSEL
#include "atomsel/vmd.hxx"
#include "atomsel/msys_keyword.hxx"
#endif

namespace desres { namespace msys { 

    IdList Atomselect(SystemPtr ptr, const std::string& txt) {
#ifndef MSYS_WITHOUT_ATOMSEL
        atomsel::vmd::StrList s;
        atomsel::PredicatePtr pred = atomsel::vmd::parse(txt,ptr,s);
        atomsel::Selection sel = atomsel::full_selection(ptr);
        pred->eval(sel);
        return sel.ids();
#else
        MSYS_FAIL("No atom selection support");
#endif
    }

}}
