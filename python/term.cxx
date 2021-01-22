#include "pymod.hxx"
#include <pybind11/stl.h>
#include "capsule.hxx"

#include <msys/term_table.hxx>
#include <msys/override.hxx>

namespace {

    Id add_term(TermTable& table, IdList const&ids, object param) {
        Id id=BadId;
        if (!param.is_none()) id = param.cast<Id>();
        return table.addTerm(ids, id);
    }

    void set_param(TermTable& table, Id term, object param) {
        Id id=BadId;
        if (!param.is_none()) id = param.cast<Id>();
        table.setParam(term, id);
    }

    handle term_prop_type(TermTable& table, Id col) {
        return from_value_type(table.termPropType(col));
    }
    Id add_term_prop( TermTable& table, String const& name, object type ) {
        return table.addTermProp(name, as_value_type(type));
    }

    /* TODO: permit col to be provided as string */
    object get_term_prop(TermTable& p, Id row, Id col) {
        return from_value_ref(p.termPropValue(row,col));
    }
    void set_term_prop(TermTable& p, Id row, Id col, object newval) {
        to_value_ref(newval, p.termPropValue(row,col));
    }

    object get_prop(TermTable& p, Id row, Id col) {
        return from_value_ref(p.propValue(row,col));
    }
    void set_prop(TermTable& p, Id row, Id col, object newval) {
        to_value_ref(newval, p.propValue(row,col));
    }

}

namespace desres { namespace msys { 

    void export_term(module m) {

        enum_<Category>(m, "Category")
            .value("none", NO_CATEGORY)
            .value("bond", BOND)
            .value("constraint", CONSTRAINT)
            .value("virtual", VIRTUAL)
            .value("polar", POLAR)
            .value("nonbonded", NONBONDED)
            .value("exclusion", EXCLUSION)
            ;

        m.def("parse_category", parse);
        m.def("print_category", print);

        class_<TermTable, TermTablePtr>(m, "TermTablePtr")
            .def("__eq__", [](TermTable* self, TermTable* other) { return self==other; })
            .def("__ne__", [](TermTable* self, TermTable* other) { return self!=other; })
            .def("__hash__", [](TermTable* g) { return size_t(g); })
            .def("destroy",     &TermTable::destroy)
            .def("system",      &TermTable::system)
            .def("atomCount",   &TermTable::atomCount)
            .def("termCount",   &TermTable::termCount)
            .def_readwrite("category", &TermTable::category)
            .def("name",        &TermTable::name)
            .def("rename",      &TermTable::rename)
            .def("terms",       &TermTable::terms)
            .def("addTerm",     add_term)
            .def("hasTerm",     &TermTable::hasTerm)
            .def("delTerm",     &TermTable::delTerm)
            .def("atoms",       &TermTable::atoms)
            .def("atom",        &TermTable::atom)
            .def("params",      &TermTable::params)
            .def("param",       &TermTable::param)
            .def("setParam",    set_param)

            /* param properties */
            .def("getProp",     get_prop)
            .def("setProp",     set_prop)

            /* term properties */
            .def("termPropCount",&TermTable::termPropCount)
            .def("termPropName", &TermTable::termPropName)
            .def("termPropIndex",&TermTable::termPropIndex)
            .def("termPropType", term_prop_type)
            .def("addTermProp",  add_term_prop)
            .def("delTermProp",  &TermTable::delTermProp)
            .def("getTermProp",  get_term_prop)
            .def("setTermProp",  set_term_prop)

            /* lookup terms based on atom */
            .def("delTermsWithAtom",    &TermTable::delTermsWithAtom)
            .def("findWithAll", &TermTable::findWithAll)
            .def("findWithAny", &TermTable::findWithAny)
            .def("findWithOnly",&TermTable::findWithOnly)
            .def("findExact"  , &TermTable::findExact)

            /* coalesce */
            .def("coalesce",    &TermTable::coalesce)

            /* overrides */
            .def("overrides",   &TermTable::overrides)

            /* misc */
            .def("resetParams", &TermTable::resetParams)
            .def("replaceWithSortedTerms", ReplaceTableWithSortedTerms)
            .def_static("asCapsule", python::termtable_as_capsule)
            .def_static("fromCapsule", python::termtable_from_capsule)
            ;
    }
}}
