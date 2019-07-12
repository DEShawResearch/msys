#include "wrap_obj.hxx"
#include "value.hxx"

using namespace boost::python;
using namespace desres::msys;

namespace {
#if PY_MAJOR_VERSION >= 3
    auto py_from_long = PyLong_FromLong;
    auto py_from_string = PyUnicode_FromString;
#else
    auto py_from_long = PyInt_FromLong;
    auto py_from_string = PyString_FromString;
#endif
}


namespace {

    object propmap_keys(VariantMap& p) {
        list s;
        for (VariantMap::value_type keyvals : p) {
            s.append(keyvals.first);
        }
        return std::move(s);
    }

    class get_visitor : public boost::static_visitor<> {
        PyObject **o;
    public:
        explicit get_visitor(PyObject** o) : o(o) {}
        void operator()(Int& i) const { *o = py_from_long(i); }
        void operator()(Float& i) const { *o = PyFloat_FromDouble(i); }
        void operator()(String& i) const { *o = py_from_string(i.data()); }
    };

    PyObject* propmap_get(VariantMap& p, std::string const& key) {
        VariantMap::iterator i=p.find(key);
        if (i==p.end()) {
            Py_INCREF(Py_None);
            return Py_None;
        }
        PyObject* o = NULL;
        boost::apply_visitor(get_visitor(&o), i->second);
        return o;
    }

    void propmap_set(VariantMap& p, std::string const& key, object type, 
                     object val) {
        ValueType t = as_value_type(type);
        Variant& v = p[key];
        switch (t) {
            case IntType: v = extract<int64_t>(val); break;
            case FloatType: v = extract<double>(val); break;
            case StringType:
            default:        v = extract<std::string>(val); break;
        }
    }
    void propmap_del(VariantMap& p, std::string const& key) {
        VariantMap::iterator i=p.find(key);
        if (i!=p.end()) p.erase(i);
    }
}

namespace desres { namespace msys { 

    void export_propmap() {
        class_<VariantMap>("VariantMap", no_init)
            .def("empty", &VariantMap::empty)
            .def("keys",  propmap_keys)
            .def("get",   propmap_get)
            .def("set",   propmap_set)
            .def("_del",  propmap_del)
            ;
    }

}}

