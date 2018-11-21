#include <boost/python.hpp>
#include "system.hxx"

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

    template <class Obj>
    unsigned long obj_hash( Obj const& obj ) { 
        return reinterpret_cast<unsigned long>(obj.get());
    }
}

namespace desres { namespace msys {

    /* TODO: make a to-python converter instead */
    object from_value_ref(const ValueRef& val);

    /* TODO: make a from-python converter instead */
    void to_value_ref(object& obj, ValueRef val);

    ValueType as_value_type(object& typeobj);
    PyObject* from_value_type(ValueType type);

    list to_python(MultiIdList const& m);
    list to_python(IdList const& m);

    IdList ids_from_python(list m);

}}

