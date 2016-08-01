#include "wrap_obj.hxx"
#include "param_table.hxx"

using namespace desres::msys;

namespace {

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

    list find_int(ParamTable& p, Id col, Int val) {
        return to_python(p.findInt(col,val));
    }
    list find_float(ParamTable& p, Id col, Float val) {
        return to_python(p.findFloat(col,val));
    }
    list find_string(ParamTable& p, Id col, String const& val) {
        return to_python(p.findString(col,val));
    }

    list wrap_params(ParamTable& p) {
        return to_python(p.params());
    }

}

namespace desres { namespace msys { 

    void export_param() {

        class_<ParamTable, ParamTablePtr>("ParamTablePtr", no_init)
            .def("__eq__",      list_eq<ParamTablePtr>)
            .def("__ne__",      list_ne<ParamTablePtr>)
            .def("__hash__",   obj_hash<ParamTablePtr>)
            .def("create",     &ParamTable::create).staticmethod("create")
            .def("paramCount", &ParamTable::paramCount)
            .def("addParam",   &ParamTable::addParam)
            .def("hasParam",   &ParamTable::hasParam)
            .def("propCount",  &ParamTable::propCount)
            .def("propName",   &ParamTable::propName)
            .def("propType",   param_prop_type)
            .def("propIndex",  &ParamTable::propIndex)
            .def("addProp",    param_add_prop)
            .def("delProp",    &ParamTable::delProp)
            .def("params",     wrap_params)
            .def("getProp",    param_get_value)
            .def("setProp",    param_set_value)
            .def("duplicate",  &ParamTable::duplicate)
            .def("refcount",   &ParamTable::refcount)
            .def("compare",    &ParamTable::compare)
            .def("findInt",    find_int)
            .def("findFloat",  find_float)
            .def("findString", find_string)
            ;
    }

}}
