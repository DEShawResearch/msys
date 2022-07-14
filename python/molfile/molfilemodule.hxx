#ifndef MOLFILE_MODULE_HXX
#define MOLFILE_MODULE_HXX

#include <pybind11/pybind11.h>
#include <msys/molfile/libmolfile_plugin.h>
#include <msys/molfile/molfile.hxx>

using namespace pybind11;

// This python module exposes the VMD plugin API through a custom 
// type called 'plugin'.  Instances of this type correspond to the
// objects contained in the module.

namespace desres { namespace molfile {


    void export_plugin(module m);
    void export_reader(module m);
    void export_frame(module m);
    void export_writer(module m);

    struct Atom_t {
        PyObject_HEAD
            molfile_atom_t atom;
        PyObject * bonds;   // dictionary from Atom_t to bond properties
        int optflags;
    };

    extern PyTypeObject AtomType;
    int initialize_atom();

    enum DataType { INT, FLOAT, DOUBLE };

    // return a NumPy array of the given shape and type.  If base is NULL, the
    // data will be copied to the numpy array and the array takes ownership.
    // If base is non-NULL, the data is not copied, and base will be the owner
    // of the data.
    PyObject **object_array(int size, PyObject **returned_result);

}}

#endif
