#include "pymod.hxx"
#include "system.hxx"
#include "version.hxx"

using namespace pybind11;

namespace desres { namespace msys {
    void export_analyze(module m);
    void export_annotated_system(module m);
    void export_chain(module m);
    void export_param(module m);
    void export_system(module m);
    void export_term(module m);
    void export_override(module m);
    void export_io(module m);
    void export_graph(module m);
    void export_inchi(module m);
    void export_spatial_hash(module m);
    void export_ff(module m);
}}

PYBIND11_MODULE(_msys, m) {
    if (MSYS_ABI_VERSION != desres::msys::abi_version()) {
        PyErr_Format(PyExc_RuntimeError,
                "This module was compiled with msys ABI version %d, but a package using msys version %d was loaded first.",
                MSYS_ABI_VERSION, desres::msys::abi_version());
        throw error_already_set();
    }
    desres::msys::export_analyze(m);
    desres::msys::export_annotated_system(m);
    desres::msys::export_system(m);
    desres::msys::export_chain(m);
    desres::msys::export_param(m);
    desres::msys::export_term(m);
    desres::msys::export_override(m);
    desres::msys::export_io(m);
    desres::msys::export_graph(m);
    desres::msys::export_inchi(m);
    desres::msys::export_spatial_hash(m);
    desres::msys::export_ff(m);
}

namespace desres { namespace msys {

    object from_value_ref(const ValueRef& val) {
        switch (val.type()) {
            case IntType: return int_(val.asInt());
            case FloatType: return float_(val.asFloat());
            case StringType: return str(val.asString());
            default:;
        }
        return none();
    }

    void to_value_ref(object obj, ValueRef val) {
        switch (val.type()) {
            case IntType:
                val.fromInt(obj.cast<Int>());
                break;
            case FloatType:
                val.fromFloat(obj.cast<Float>());
                break;
            case StringType:
                val.fromString(obj.cast<String>());
        }
    }

    ValueType as_value_type(handle typeobj) {
        auto ptr = reinterpret_cast<char *>(typeobj.ptr());
        if (ptr == (char *)&PyFloat_Type) return FloatType;
        if (ptr == (char *)&PyLong_Type) return IntType;
        if (ptr == (char *)&PyUnicode_Type) return StringType;
        PyErr_SetString(PyExc_ValueError, "Expected int, float or str as type argument");
        throw error_already_set();
    }

    handle from_value_type(ValueType type) {
        switch (type) {
            case IntType: return (PyObject *)&PyLong_Type;
            case FloatType: return (PyObject *)&PyFloat_Type;
            case StringType: return (PyObject *)&PyUnicode_Type;
            default:;
        }
        return none();
    }

}}

