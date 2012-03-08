#include "atomsel.hxx"
#include "atomsel/vmd.hxx"
#include "atomsel/msys_keyword.hxx"

namespace desres { namespace msys { 

    IdList Atomselect(SystemPtr ptr, const std::string& txt) {
        atomsel::vmd::StrList s;
        atomsel::PredicatePtr pred = atomsel::vmd::parse(txt,ptr,s);
        atomsel::Selection sel = atomsel::full_selection(ptr);
        pred->eval(sel);
        return sel.ids();
    }

}}
