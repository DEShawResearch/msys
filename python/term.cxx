#include "wrap_obj.hxx"
#include "term_table.hxx"

using namespace desres::msys;

namespace {

    Id add_term(TermTable& table, object atoms, object param) {
        Id id=BadId;
        if (param.ptr()!=Py_None) id=extract<Id>(param);
        list L(atoms);
        IdList ids(len(L));
        for (unsigned i=0; i<ids.size(); i++) ids[i]=extract<Id>(L[i]);
        return table.addTerm(ids, id);
    }

    void set_param(TermTable& table, Id term, object param) {
        Id id=BadId;
        if (param.ptr()!=Py_None) id=extract<Id>(param);
        table.setParam(term,id);
    }

    void set_paramB(TermTable& table, Id term, object param) {
        Id id=BadId;
        if (param.ptr()!=Py_None) id=extract<Id>(param);
        table.setParamB(term,id);
    }
    
    PyObject* term_prop_type(TermTable& table, Id col) {
        return from_value_type(table.termPropType(col));
    }
    Id add_term_prop( TermTable& table, String const& name, object type ) {
        return table.addTermProp(name, as_value_type(type));
    }
    object get_term_prop(TermTable& p, Id row, Id col) {
        return from_value_ref(p.termPropValue(row,col));
    }
    void set_term_prop(TermTable& p, Id row, Id col, object newval) {
        to_value_ref(newval, p.termPropValue(row,col));
    }
    PyObject* prop_type(TermTable& table, Id col) {
        return from_value_type(table.propType(col));
    }
    Id add_prop( TermTable& table, String const& name, object type ) {
        return table.addProp(name, as_value_type(type));
    }
    object get_prop(TermTable& p, Id row, Id col) {
        return from_value_ref(p.propValue(row,col));
    }
    void set_prop(TermTable& p, Id row, Id col, object newval) {
        to_value_ref(newval, p.propValue(row,col));
    }

}

namespace desres { namespace msys { 

    void export_term() {

        class_<TermTable, TermTablePtr>("TermTablePtr", no_init)
            .def("__eq__",      list_eq<TermTablePtr>)
            .def("__ne__",      list_ne<TermTablePtr>)
            .def("system",      &TermTable::system)
            .def("atomCount",   &TermTable::atomCount)
            .def("termCount",   &TermTable::termCount)
            .def_readwrite("category", &TermTable::category)
            .def("terms",       &TermTable::terms)
            .def("addTerm",     add_term)
            .def("hasTerm",     &TermTable::hasTerm)
            .def("delTerm",     &TermTable::delTerm)
            .def("atoms",       &TermTable::atoms)
            .def("paramTable",  &TermTable::paramTable)
            .def("param",       &TermTable::param)
            .def("setParam",    set_param)
            .def("paramB",      &TermTable::paramB)
            .def("setParamB",   set_paramB)
            .def("alchemical",  &TermTable::alchemical)

            /* term properties */
            .def("termPropCount",&TermTable::termPropCount)
            .def("termPropName", &TermTable::termPropName)
            .def("termPropIndex",&TermTable::termPropIndex)
            .def("termPropType", term_prop_type)
            .def("addTermProp",  add_term_prop)
            .def("getTermProp",  get_term_prop)
            .def("setTermProp",  set_term_prop)

            /* parameter properties */
            .def("propCount",   &TermTable::propCount)
            .def("propName",    &TermTable::propName)
            .def("propIndex",   &TermTable::propIndex)
            .def("propType",    prop_type)
            .def("addProp",  add_prop)
            .def("getProp",  get_prop)
            .def("setProp",  set_prop)
            ;
    }

}}
