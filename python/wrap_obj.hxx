#include <boost/python.hpp>
#include "system.hxx"
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

using namespace boost::python;

namespace {

    using namespace desres::msys;

    typedef with_custodian_and_ward_postcall<0,1> return_obj;
    typedef return_internal_reference<> return_ptr;
    typedef return_value_policy<copy_const_reference> return_const;

    template <class List>
    bool list_eq(const List& self, const List& other) { return self==other; }

    template <class List>
    bool list_ne(const List& self, const List& other) { return self!=other; }

    template <class List>
    class_<List> declare_list( const char * name ) {
        return class_<List>(name)
            .def(vector_indexing_suite<List>())
            /* for some reason, operator== returns false even when all
             * the component elements compare equal.  What's that about? */
            .def("__eq__", list_eq<List>)
            .def("__ne__", list_ne<List>)
            ;
    
    }

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

}

