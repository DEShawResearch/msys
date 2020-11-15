#include <pybind11/pybind11.h>
#include <msys/value.hxx>

using namespace desres::msys;
using namespace pybind11;

namespace desres { namespace msys {


    object from_value_ref(const ValueRef& val);
    void to_value_ref(object obj, ValueRef val);
    ValueType as_value_type(handle typeobj);
    handle from_value_type(ValueType type);

}}

