#include "wrap_obj.hxx"

namespace {
    list ct_keys(component_t& ct) {
        list L;
        BOOST_FOREACH(String const& key, ct.keys()) {
            L.append(key);
        }
        return L;
    }
    object ct_get(component_t& ct, String const& key) {
        if (!ct.has(key)) {
            PyErr_Format(PyExc_KeyError, "No property '%s' in ct", key.c_str());
            throw_error_already_set();
        }
        return from_value_ref(ct.value(key));
    }

    void ct_add(component_t& ct, String const& key, object obj) {
        ct.add(key,as_value_type(obj));
    }

    void ct_set(component_t& ct, String const& key, object obj) {
        to_value_ref(obj, ct.value(key));
    }
    
    PyObject* ct_type(component_t& ct, String const& key) {
        if (!ct.has(key)) {
            Py_INCREF(Py_None);
            return Py_None;
        }
        return from_value_type(ct.type(key));
    }
}

namespace desres { namespace msys { 

    void export_chain() {

        class_<chain_t>("chain_t", no_init)
            .def_readonly("ct", &chain_t::ct)
            .def_readwrite("name", &chain_t::name)
            .def_readwrite("segid", &chain_t::segid)
            ;

        class_<component_t>("component_t", no_init)
            .add_property("name", &component_t::name, &component_t::setName)
            .def("keys",ct_keys)
            .def("get", ct_get)
            .def("add", ct_add)
            .def("set", ct_set)
            .def("type", ct_type)
            .def("remove", &component_t::del)
            ;
    }
}}

