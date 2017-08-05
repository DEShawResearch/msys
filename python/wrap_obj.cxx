#include "wrap_obj.hxx"

using namespace boost::python;

namespace {
#if PY_MAJOR_VERSION >= 3
    auto py_int_type = &PyLong_Type;
    auto py_str_type = &PyUnicode_Type;
#else
    auto py_int_type = &PyInt_Type;
    auto py_str_type = &PyString_Type;
#endif
}

namespace desres { namespace msys {

    /* TODO: make a to-python converter instead */
    object from_value_ref(const ValueRef& val) {
        if (val.type()==IntType) return object(val.asInt());
        if (val.type()==FloatType) return object(val.asFloat());
        if (val.type()==StringType) return object(val.asString());
        return object();
    }

    /* TODO: make a from-python converter instead */
    void to_value_ref(object& obj, ValueRef val) {
        if (val.type()==IntType) val.fromInt(extract<Int>(obj));
        if (val.type()==FloatType) val.fromFloat(extract<Float>(obj));
        if (val.type()==StringType) val.fromString(extract<String>(obj));
    }

    ValueType as_value_type(object& typeobj) {
        ValueType type = IntType;   /* silence compiler warning */
        if ((char *)typeobj.ptr()==(char *)&PyFloat_Type) {
            type=FloatType;
        } else if ((char *)typeobj.ptr()==(char *)py_int_type) {
            type=IntType;
        } else if ((char *)typeobj.ptr()==(char *)py_str_type) {
            type=StringType;
        } else {
            PyErr_Format(PyExc_ValueError,
                    "Expected int, float or string as type argument");
            throw_error_already_set();
        }
        return type;
    }

    PyObject* from_value_type(ValueType type) {
        PyTypeObject* typeobj=NULL;
        switch (type) {
            case IntType: typeobj=py_int_type; break;
            case FloatType: typeobj=&PyFloat_Type; break;
            case StringType: typeobj=py_str_type; break;
            default: throw std::runtime_error("unrecognized type");
        }
        PyObject* obj = (PyObject*)typeobj;
        Py_INCREF(obj);
        return obj;
    }

    list to_python(IdList const& mm) {
        list result;
        for (auto id : mm) result.append(id);
        return result;
    }

    list to_python(MultiIdList const& m) {
        list result;
        for (auto const& mm : m)  result.append(to_python(mm));
        return result;
    }

    IdList ids_from_python(list m) {
        IdList ids;
        ids.reserve(len(m));
        for (Id i=0, n=len(m); i<n; i++) {
            ids.push_back(extract<Id>(m[i]));
        }
        return ids;
    }

}}

