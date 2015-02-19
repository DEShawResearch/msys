#ifndef MOLFILE_MODULE_HXX
#define MOLFILE_MODULE_HXX

#include <Python.h>
#include "molfile/libmolfile_plugin.h"
#include "molfile/molfile.hxx"

// This python module exposes the VMD plugin API through a custom 
// type called 'plugin'.  Instances of this type correspond to the
// objects contained in the module.

namespace desres { namespace molfile {

    void export_plugin();
    void export_reader();
    void export_frame();
    void export_writer();
    void export_dtrreader();

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
    PyObject *backed_vector( int nd, Py_ssize_t *dims, DataType type, void *data, PyObject *base );
    PyObject **object_array(int size, PyObject **returned_result);
    void * array_data( PyObject * arr );

}}

#endif
