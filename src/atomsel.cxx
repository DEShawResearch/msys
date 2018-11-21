#include "atomsel.hxx"
#include "atomsel/token.hxx"

namespace desres { namespace msys { 

    IdList Atomselect(SystemPtr ptr, const std::string& txt) {
        return Atomselect(ptr, txt, nullptr, nullptr);
    }

    IdList Atomselect(SystemPtr ptr, const std::string& txt,
                      const float* pos, const double* cell) {

        atomsel::Query q;
        q.mol = ptr.get();
        q.pos = pos;
        q.cell = cell;
        q.parse(txt);
        auto s = atomsel::full_selection(q.mol);
        q.pred->eval(s);
        return s.ids();
    }

}}
