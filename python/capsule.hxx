#ifndef desres_msys_python_capsule_hxx
#define desres_msys_python_capsule_hxx

#include <Python.h>
#include <msys/system.hxx>

namespace desres { namespace msys { namespace python {

    // Functions for converting boost::python wrapped msys classes
    // to and from Python capsule objects.

    // tag for msys.System.  Don't use directly.
    static const char* _system_capsule_name = "_msys.SystemPtr";

    // construct a Python capsule from a SystemPtr.  The lifetime of the
    // system must be at least as long as the capsule; if the system 
    // refcount drops to zero the pointer returned by the capsule will
    // be invalid! 
    static inline
    PyObject* system_as_capsule(SystemPtr mol) {
        return PyCapsule_New(mol.get(), _system_capsule_name, nullptr);
    }

    // Extract a SystemPtr from a Python capsule assumed to have been
    // created from system_as_capsule().  Returns null pointer if the
    // capsule has the wrong tag name.  Throws if the pointer held by
    // the capsule was not created as a shared pointer.  Undefined behavior
    // if the capsule holds the wrong kind of pointer, or the pointer
    // has gone out of scope.
    static inline
    SystemPtr system_from_capsule(PyObject* obj) {
        void* ptr = PyCapsule_GetPointer(obj, _system_capsule_name);
        if (ptr == nullptr) return SystemPtr();
        return reinterpret_cast<System*>(ptr)->shared_from_this();
    }

    /* Example using pybind11:
#include <msys/python/capsule.hxx>
using namespace desres::msys::python;

object system_inout(object obj) {
    auto mol = system_from_capsule(obj.ptr());
    if (!mol) throw error_already_set();
    std::cout << "got mol " << mol->name
              << " with " << mol->atomCount() << " atoms.\n";
    return reinterpret_steal<object>(system_as_capsule(mol));
}
    */

}}}

#endif
