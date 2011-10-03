#include "wrap_obj.hxx"
#include "param_table.hxx"

using namespace desres::msys;

namespace {

    template <typename T>
    bool ptr_eq(boost::shared_ptr<T> p, boost::shared_ptr<T> q) {
        return p==q;
    }

    Id param_add_prop(ParamTable& p, const std::string& name, object typeobj) {
        return p.addProp(name, as_value_type(typeobj));
    }
    object param_get_value(ParamTable& p, Id row, Id col) {
        return from_value_ref(p.value(row,col));
    }
    void param_set_value(ParamTable& p, Id row, Id col, object newval) {
        to_value_ref(newval, p.value(row,col));
    }

    PyObject* param_prop_type(ParamTable& p, Id col) {
        return from_value_type(p.propType(col));
    }

}

namespace desres { namespace msys { 

    void export_param() {

        class_<ParamTable, ParamTablePtr>("ParamTablePtr", no_init)
            .def("__eq__",      list_eq<ParamTablePtr>)
            .def("__ne__",      list_ne<ParamTablePtr>)
            .def("create",     &ParamTable::create).staticmethod("create")
            .def("paramCount", &ParamTable::paramCount)
            .def("addParam",   &ParamTable::addParam)
            .def("hasParam",   &ParamTable::hasParam)
            .def("propCount",  &ParamTable::propCount)
            .def("propName",   &ParamTable::propName)
            .def("propType",   param_prop_type)
            .def("propIndex",  &ParamTable::propIndex)
            .def("addProp",    param_add_prop)
            .def("ids",        &ParamTable::params)
            .def("getProp",    param_get_value)
            .def("setProp",    param_set_value)
            ;
    }

}}
