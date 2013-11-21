#include "wrap_obj.hxx"

using namespace boost::python;

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
        } else if ((char *)typeobj.ptr()==(char *)&PyInt_Type) {
            type=IntType;
        } else if ((char *)typeobj.ptr()==(char *)&PyString_Type) {
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
            case IntType: typeobj=&PyInt_Type; break;
            case FloatType: typeobj=&PyFloat_Type; break;
            case StringType: typeobj=&PyString_Type; break;
            default: throw std::runtime_error("unrecognized type");
        }
        PyObject* obj = (PyObject*)typeobj;
        Py_INCREF(obj);
        return obj;
    }

}}

