#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <msys/override.hxx>

using namespace desres::msys;
using namespace pybind11;

namespace desres { namespace msys { 

    void export_override(module m) {

        class_<OverrideTable, OverrideTablePtr>(m, "OverrideTablePtr")
            .def("__eq__", [](OverrideTable* self, OverrideTable* other) { return self==other; })
            .def("__ne__", [](OverrideTable* self, OverrideTable* other) { return self!=other; })
            .def("__hash__", [](OverrideTable* g) { return size_t(g); })
            .def("target",      &OverrideTable::target)
            .def("params",      &OverrideTable::params)
            .def("get", [](OverrideTable& o, Id p1, Id p2) { return o.get({p1,p2}); })
            .def("del_",[](OverrideTable& o, Id p1, Id p2) { o.del({p1,p2}); })
            .def("set",[](OverrideTable& o, Id p1, Id p2, Id p) { o.set({p1,p2}, p); })
            .def("count", &OverrideTable::count)
            .def("list", &OverrideTable::list)
            ;
    }

}}
