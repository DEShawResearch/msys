#include "pymod.hxx"
#include <pybind11/stl.h>
#include <msys/param_table.hxx>
#include "capsule.hxx"

namespace desres { namespace msys { 

    void export_param(module m) {

        class_<ParamTable, ParamTablePtr>(m, "ParamTablePtr")
            .def("__eq__", [](ParamTable* self, ParamTable* other) { return self==other; })
            .def("__ne__", [](ParamTable* self, ParamTable* other) { return self!=other; })
            .def("__hash__", [](ParamTable* g) { return size_t(g); })
            .def_static("create", &ParamTable::create)
            .def("paramCount", &ParamTable::paramCount)
            .def("addParam",   &ParamTable::addParam)
            .def("hasParam",   &ParamTable::hasParam)
            .def("propCount",  &ParamTable::propCount)
            .def("propName",   &ParamTable::propName)
            .def("propType",   [](ParamTable& p, Id id) { return from_value_type(p.propType(id)); })
            .def("propIndex",  &ParamTable::propIndex)
            .def("addProp",    [](ParamTable& t, String const& name, object type) { t.addProp(name, as_value_type(type)); })
            .def("delProp",    &ParamTable::delProp)
            .def("params",     &ParamTable::params)
            .def("getProp",    [](ParamTable& t, Id row, Id col) { return from_value_ref(t.value(row,col)); })
            .def("setProp",    [](ParamTable& t, Id row, Id col, object val) { to_value_ref(val, t.value(row,col)); })
            .def("duplicate",  &ParamTable::duplicate)
            .def("refcount",   &ParamTable::refcount)
            .def("compare",    &ParamTable::compare)
            .def("findInt",    &ParamTable::findInt)
            .def("findFloat",  &ParamTable::findFloat)
            .def("findString", &ParamTable::findString)
            .def_static("asCapsule", [](ParamTablePtr p) { return handle(python::paramtable_as_capsule(p)); })
            .def_static("fromCapsule", [](handle h) { return python::paramtable_from_capsule(h.ptr()); })
            .def("valuesForColumn", [](ParamTablePtr p, Id col) {
                switch (p->propType(col)) {
                case IntType: return cast(p->valuesForColumn<Int>(col));
                case FloatType: return cast(p->valuesForColumn<Float>(col));
                default:
                case StringType: return cast(p->valuesForColumn<String>(col));
                }; });
    }

}}
