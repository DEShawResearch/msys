#include "molfilemodule.hxx"

using namespace desres::molfile;

namespace {
    auto py_as_long = PyLong_AsLong;
    auto py_from_long = PyLong_FromLong;
    auto py_as_string = PyUnicode_AsUTF8;
    auto py_from_string = PyUnicode_FromString;
    auto py_from_format = PyUnicode_FromFormat;
}

#define GETSET_INT(ATTR, FLAG) \
    PyObject *get_##ATTR( PyObject *pySelf, void *closure) { \
        Atom_t *self = reinterpret_cast<Atom_t*>(pySelf); \
        return py_from_long(self->atom.ATTR); \
    } \
\
int set_##ATTR( PyObject *pySelf, PyObject *rval, void *closure) { \
    Atom_t *self = reinterpret_cast<Atom_t*>(pySelf); \
    int val=py_as_long(rval); \
    if (PyErr_Occurred()) return -1; \
    self->atom.ATTR = val; \
    self->optflags |= FLAG; \
    return 0; \
}

#define GETSET_FLT(ATTR, FLAG) \
    PyObject *get_##ATTR( PyObject *pySelf, void *closure) { \
        Atom_t *self = reinterpret_cast<Atom_t*>(pySelf); \
        return PyFloat_FromDouble(self->atom.ATTR); \
    } \
\
int set_##ATTR( PyObject *pySelf, PyObject *rval, void *closure) { \
    Atom_t *self = reinterpret_cast<Atom_t*>(pySelf); \
    double val=PyFloat_AsDouble(rval); \
    if (PyErr_Occurred()) return -1; \
    self->atom.ATTR = val; \
    self->optflags |= FLAG; \
    return 0; \
}

#define GETSET_STR(ATTR, FLAG) \
    PyObject *get_##ATTR( PyObject *pySelf, void *closure) { \
        Atom_t *self = reinterpret_cast<Atom_t*>(pySelf); \
        return py_from_string(self->atom.ATTR); \
    } \
\
int set_##ATTR( PyObject *pySelf, PyObject *rval, void *closure) { \
    Atom_t *self = reinterpret_cast<Atom_t*>(pySelf); \
    const char *val=py_as_string(rval); \
    if (PyErr_Occurred()) return -1; \
    strncpy(self->atom.ATTR, val, sizeof(self->atom.ATTR)); \
    self->atom.ATTR[sizeof(self->atom.ATTR)-1]='\0'; \
    self->optflags |= FLAG; \
    return 0; \
}

namespace {

    int atom_init(PyObject *pySelf, PyObject *args, PyObject *kwds) {

        // no positional arguments
        Py_ssize_t nargs=PySequence_Size(args);
        if (nargs!=0) {
            PyErr_Format(PyExc_TypeError, "Atom() takes no positional arguments (%ld given", nargs);
            return -1;
        }
        Atom_t *self = reinterpret_cast<Atom_t*>(pySelf);
        memset(&self->atom, 0, sizeof(self->atom));
        self->atom.chain[0]=' ';
        self->bonds = PyDict_New();

        // take all keyword arguments and set as attributes
        if (kwds) {
            PyObject* key, *val;
            Py_ssize_t pos = 0;
            while (PyDict_Next(kwds, &pos, &key, &val)) {
                if (PyObject_SetAttr(pySelf, key, val)) {
                    return -1;
                }
            }
        }
        return 0;
    }

    void atom_dealloc(PyObject *pySelf) {
        PyObject_GC_UnTrack(pySelf);
        Atom_t *self = reinterpret_cast<Atom_t*>(pySelf);
        Py_XDECREF(self->bonds); // might have been decref'd by Py_CLEAR
        PyObject_GC_Del(pySelf);
    }

    GETSET_STR(name,0)
    GETSET_STR(type,0)
    GETSET_STR(resname,0)
    GETSET_INT(resid,0)
    GETSET_STR(segid,0)
    GETSET_STR(chain,0)

    GETSET_STR(altloc, MOLFILE_ALTLOC)
    GETSET_STR(insertion, MOLFILE_INSERTION)

    GETSET_FLT(occupancy, MOLFILE_OCCUPANCY)
    GETSET_FLT(bfactor, MOLFILE_BFACTOR)
    GETSET_FLT(mass, MOLFILE_MASS)
    GETSET_FLT(charge, MOLFILE_CHARGE)
    GETSET_FLT(radius, MOLFILE_RADIUS)

    GETSET_INT(atomicnumber, MOLFILE_ATOMICNUMBER)

    PyObject *atom_repr(PyObject *pySelf) {
        Atom_t *self = reinterpret_cast<Atom_t*>(pySelf);
        return py_from_format("<Atom '%s'>", self->atom.name);
    }

    PyObject *get_bonds(PyObject *pySelf) {
        Atom_t *self = reinterpret_cast<Atom_t*>(pySelf);
        return PyFrozenSet_New(self->bonds);
    }

    PyGetSetDef atom_getset[] = {
        { (char *)"bonds",(getter)get_bonds, NULL, (char *)"bonded atoms", NULL },
        { (char *)"name", (getter)get_name, set_name, (char *)"atom name", NULL },
        { (char *)"type", (getter)get_type, set_type, (char *)"atom type", NULL },
        { (char *)"resname", (getter)get_resname, set_resname, (char *)"residue name", NULL },
        { (char *)"resid", (getter)get_resid, set_resid, (char *)"residue id", NULL },
        { (char *)"segid", (getter)get_segid, set_segid, (char *)"segment name", NULL },
        { (char *)"chain", (getter)get_chain, set_chain, (char *)"chain name", NULL },

        { (char *)"altloc", (getter)get_altloc, set_altloc, (char *)"segment name", NULL },
        { (char *)"insertion", (getter)get_insertion, set_insertion, (char *)"segment name", NULL },

        { (char *)"occupancy", (getter)get_occupancy, set_occupancy, (char *)"occupancy", NULL },
        { (char *)"bfactor", (getter)get_bfactor, set_bfactor, (char *)"temperature factor", NULL },
        { (char *)"mass", (getter)get_mass, set_mass, (char *)"mass in amu", NULL },
        { (char *)"charge", (getter)get_charge, set_charge, (char *)"charge", NULL },
        { (char *)"radius", (getter)get_radius, set_radius, (char *)"vdw radius", NULL },

        { (char *)"anum", (getter)get_atomicnumber, set_atomicnumber, (char *)"atomic number", NULL },
        { NULL },
    };

    const char *addbond_doc =
        "addbond(atom) -- add atom to self.bonds and self to atom.bonds\n";
    PyObject *addbond(PyObject *pySelf, PyObject *args, PyObject *kwds) {
        static char *kwlist[] = { (char *)"atom", (char *)"order", 0 };
        PyObject *pyAtom=NULL;
        float order = 1;
        if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|f", kwlist,
                    &AtomType, &pyAtom, &order))
            return NULL;
        Atom_t *self = reinterpret_cast<Atom_t*>(pySelf);
        Atom_t *atom = reinterpret_cast<Atom_t*>(pyAtom);
        if (atom == self) {
            PyErr_SetString(PyExc_ValueError, "Cannot add bond to self");
            return NULL;
        }
        PyObject * orderobj = PyFloat_FromDouble(order);
        if (PyDict_SetItem(self->bonds, pyAtom, orderobj) ||
                PyDict_SetItem(atom->bonds, pySelf, orderobj)) {
            Py_DECREF(orderobj);
            return NULL;
        }
        Py_DECREF(orderobj);
        Py_INCREF(Py_None);
        return Py_None;
    }

    const char *delbond_doc =
        "delbond(atom) -- remove atom from self.bonds and self from atom.bonds\n";
    PyObject *delbond(PyObject *pySelf, PyObject *args, PyObject *kwds) {
        static char *kwlist[] = { (char *)"atom", 0 };
        PyObject *pyAtom=NULL;
        if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!", kwlist,
                    &AtomType, &pyAtom))
            return NULL;
        Atom_t *self = reinterpret_cast<Atom_t*>(pySelf);
        Atom_t *atom = reinterpret_cast<Atom_t*>(pyAtom);
        if (PyDict_DelItem(self->bonds, pyAtom)<0 && 
                !PyErr_ExceptionMatches(PyExc_KeyError)) 
            return NULL;
        if (PyDict_DelItem(atom->bonds, pySelf)<0 && 
                !PyErr_ExceptionMatches(PyExc_KeyError))
            return NULL;
        Py_INCREF(Py_None);
        return Py_None;
    }

    const char * getbondorder_doc =
        "getbondorder(atom) -- bond order for bond with given atom.\n";
    PyObject * getbondorder(PyObject * pySelf, PyObject * args) {
        PyObject * pyAtom=NULL;
        if (!PyArg_ParseTuple(args, "O!", &AtomType, &pyAtom))
            return NULL;
        Atom_t *self = reinterpret_cast<Atom_t*>(pySelf);
        PyObject * val = PyDict_GetItem(self->bonds, pyAtom);
        if (!val) {
            PyErr_Format(PyExc_ValueError, "no bond to atom");
            return NULL;
        }
        Py_INCREF(val);
        return val;
    }

    const char * setbondorder_doc =
        "setbondorder(atom,val) -- set bond order for bond with given atom.\n";
    PyObject * setbondorder(PyObject * pySelf, PyObject * args) {
        PyObject * pyAtom=NULL;
        float value=1;
        if (!PyArg_ParseTuple(args, "O!f", &AtomType, &pyAtom, &value))
            return NULL;
        Atom_t *self = reinterpret_cast<Atom_t*>(pySelf);
        Atom_t *atom = reinterpret_cast<Atom_t*>(pyAtom);
        PyObject * val = PyDict_GetItem(self->bonds, pyAtom);
        if (!val) {
            PyErr_Format(PyExc_ValueError, "no bond to atom");
            return NULL;
        }
        if (PyDict_SetItem(self->bonds, pyAtom, PyFloat_FromDouble(value))) {
            return NULL;
        }
        if (PyDict_SetItem(atom->bonds, pySelf, PyFloat_FromDouble(value))) {
            return NULL;
        }
        Py_INCREF(Py_None);
        return Py_None;
    }

    PyMethodDef atom_methods[] = {
        { "addbond", 
            (PyCFunction)addbond, 
            METH_VARARGS|METH_KEYWORDS,
            addbond_doc },
        { "delbond", 
            (PyCFunction)delbond, 
            METH_VARARGS|METH_KEYWORDS,
            delbond_doc },
        { "getbondorder",
            (PyCFunction)getbondorder, 
            METH_VARARGS,
            getbondorder_doc },
        { "setbondorder",
            (PyCFunction)setbondorder, 
            METH_VARARGS,
            setbondorder_doc },
        { 0,0 }
    };

    int atom_traverse(PyObject *pySelf, visitproc visit, void *arg) {
        Atom_t *self = reinterpret_cast<Atom_t*>(pySelf);
        Py_VISIT(self->bonds);
        return 0;
    }
    int atom_clear(PyObject *pySelf) {
        Atom_t *self = reinterpret_cast<Atom_t*>(pySelf);
        Py_CLEAR(self->bonds);
        return 0;
    }
}

PyTypeObject desres::molfile::AtomType;

int desres::molfile::initialize_atom() {

    AtomType.tp_name = "Atom";
    AtomType.tp_basicsize = sizeof(Atom_t);
    AtomType.tp_alloc = PyType_GenericAlloc;
    AtomType.tp_new = PyType_GenericNew;
    AtomType.tp_init = atom_init;
    AtomType.tp_dealloc = atom_dealloc;
    AtomType.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE
        | Py_TPFLAGS_HAVE_GC
        ;
    AtomType.tp_getset = atom_getset;
    AtomType.tp_methods = atom_methods;
    AtomType.tp_repr = atom_repr;
    AtomType.tp_traverse = atom_traverse;
    AtomType.tp_clear = atom_clear;

    return (PyType_Ready(&AtomType));
}

