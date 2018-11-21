#include "wrap_obj.hxx"
#include "override.hxx"

using namespace desres::msys;

namespace {
    Id get_(OverrideTable const& o, Id p1, Id p2) {
        return o.get(IdPair(p1,p2));
    }
    void del_(OverrideTable& o, Id p1, Id p2) {
        o.del(IdPair(p1,p2));
    }
    void set_(OverrideTable& o, Id p1, Id p2, Id p) {
        o.set(IdPair(p1,p2), p);
    }
    list list_(OverrideTable const& o) {
        list L;
        std::vector<IdPair> pairs(o.list());
        for (unsigned i=0; i<pairs.size(); i++) {
            list elem;
            elem.append(object(pairs[i].first));
            elem.append(object(pairs[i].second));
            L.append(elem);
        }
        return L;
    }

}

namespace desres { namespace msys { 

    void export_override() {

        class_<OverrideTable, OverrideTablePtr>("OverrideTablePtr", no_init)
            .def("__eq__",      list_eq<OverrideTablePtr>)
            .def("__ne__",      list_ne<OverrideTablePtr>)
            .def("__hash__",   obj_hash<OverrideTablePtr>)
            .def("target",      &OverrideTable::target)
            .def("params",      &OverrideTable::params)
            .def("get",         get_)
            .def("del_",        del_)
            .def("set",         set_)
            .def("count",       &OverrideTable::count)
            .def("list",        list_)
            ;
    }

}}
