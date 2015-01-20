#include "molfilemodule.hxx"
#include "findframe.hxx"
#include <vector>
#include <stdexcept>
#include <numpy/arrayobject.h>

#include <boost/python.hpp>

using namespace desres::molfile;
typedef ssize_t (Reader::*findfunc)(double T) const;

using namespace boost::python;

namespace {

    object reader_atoms(const Reader& self) {
        // build the atom list
        const std::vector<atom_t> &atoms = self.atoms();
        if (!atoms.size()) {
            PyErr_SetString(PyExc_AttributeError, "This plugin can't read atoms");
            throw error_already_set();
        }
        const std::vector<bond_t> &bonds = self.bonds();

        PyObject *result;
        PyObject **data = object_array(atoms.size(), &result);
        if (!data) throw error_already_set();
        for (unsigned i=0; i<atoms.size(); i++) {
            Atom_t *atom = PyObject_GC_New(Atom_t, &AtomType);
            if (!atom) {
                Py_DECREF(result);
                throw error_already_set();
            }
            atom->atom = atoms[i];
            atom->optflags = self.optflags();
            atom->bonds = PyDict_New();
            data[i] = reinterpret_cast<PyObject *>(atom);
            PyObject_GC_Track(data[i]);
        }
        // create bonds
        for (unsigned i=0; i<bonds.size(); i++) {
            int ai = bonds[i].from;
            int aj = bonds[i].to;
            float o = bonds[i].order;
            Atom_t *iatom = reinterpret_cast<Atom_t*>(data[ai]);
            Atom_t *jatom = reinterpret_cast<Atom_t*>(data[aj]);
            PyObject * orderobj = PyFloat_FromDouble(o);
            PyDict_SetItem( iatom->bonds, data[aj], orderobj );
            PyDict_SetItem( jatom->bonds, data[ai], orderobj );
            Py_DECREF(orderobj);
        }
        return object(handle<>(result));
    }

    object reader_topology(const Reader& self) {
        const std::vector<bond_t> &bonds = self.bonds();
        int n = self.natoms();
        PyObject * result = NULL;
        PyObject ** data = object_array(n, &result);
        if (data) {
            for (int i=0; i<n; i++) data[i] = PyTuple_New(0);
            for (unsigned i=0; i<bonds.size(); i++) {
                int from = bonds[i].from;
                int to   = bonds[i].to;
                int oldsize = PyTuple_GET_SIZE(data[from]);
                _PyTuple_Resize(&data[from], oldsize+1);
                PyTuple_SET_ITEM(data[from], oldsize, PyInt_FromLong(to));
            }
        }
        return object(handle<>(result));
    }

    object reader_bondorders(const Reader& self) {
        const std::vector<bond_t> &bonds = self.bonds();
        int n = self.natoms();
        PyObject * result = NULL;
        PyObject ** data = object_array(n, &result);
        if (data) {
            for (int i=0; i<n; i++) data[i] = PyTuple_New(0);
            for (unsigned i=0; i<bonds.size(); i++) {
                int from = bonds[i].from;
                float order = bonds[i].order;
                int oldsize = PyTuple_GET_SIZE(data[from]);
                _PyTuple_Resize(&data[from], oldsize+1);
                PyTuple_SET_ITEM(data[from], oldsize, PyFloat_FromDouble(order));
            }
        }
        return object(handle<>(result));
    }


    object reader_times(const Reader& self) {
        Py_ssize_t n = self.nframes();
        if (n<0) {
            PyErr_Format(PyExc_ValueError, 
                    "No times available for this plugin type");
            throw error_already_set();
        }
        Py_ssize_t dims[1] = { n };
        PyObject * arr = backed_vector( 1, dims, DOUBLE, NULL, NULL );
        if (self.read_times(0,n,(double *)array_data(arr))!=n) {
            Py_DECREF(arr);
            arr=NULL;
            PyErr_Format(PyExc_RuntimeError, "Error reading times");
            throw error_already_set();
        }
        return object(handle<>(arr));
    }

    template <findfunc f>
    Frame * wrap( const Reader& self, double time ) {
        Py_ssize_t index=(self.*f)(time);
        if (index<0) return NULL;
        return self.frame(index);
    }

    int reader_ngrids(Reader const& r) {
        return r.grids().size();
    }

    object reader_grid_meta(Reader const& r, int n) {
        grid_t const& g = r.grids().at(n);
        dict d;
        list origin, xaxis, yaxis, zaxis, size;
        for (int i=0; i<3; i++) {
            origin.append(g.origin[i]);
            xaxis.append(g.xaxis[i]);
            yaxis.append(g.yaxis[i]);
            zaxis.append(g.zaxis[i]);
            size.append((&g.xsize)[i]);
        }

        d["name"] = g.dataname;
        d["origin"] = origin;
        d["xaxis"] = xaxis;
        d["yaxis"] = yaxis;
        d["zaxis"] = zaxis;
        d["dims"] = size;
        return d;
    }

    void reader_grid_data(Reader const& r, int n, PyObject* arr) {
        r.read_grid(n, (float *)PyArray_DATA(arr));
    }

    Frame* reader_next(Reader& r) {
        Frame* f;
        Py_BEGIN_ALLOW_THREADS
        f = r.next();
        Py_END_ALLOW_THREADS
        return f;
    }
}

void desres::molfile::export_reader() {

    class_<Reader>("Reader", no_init)
        .add_property("natoms", &Reader::natoms)
        .add_property("nframes",&Reader::nframes)
        .add_property("ngrids", reader_ngrids)
        .add_property("has_velocities", &Reader::has_velocities)
        .add_property("atoms", reader_atoms)
        .add_property("topology", reader_topology)
        .add_property("bondorders", reader_bondorders)
        .add_property("times", reader_times)
        .def("reopen", &Reader::reopen,
                return_value_policy<manage_new_object>())
        .def("frame", &Reader::frame, 
                return_value_policy<manage_new_object>())
        .def("next", reader_next,
                return_value_policy<manage_new_object>())
        .def("skip", &Reader::skip)
        .def("at_time_near", &wrap<&Reader::at_time_near>,
                arg("time"),
                return_value_policy<manage_new_object>())
        .def("at_time_gt", &wrap<&Reader::at_time_gt>,
                arg("time"),
                return_value_policy<manage_new_object>())
        .def("at_time_ge", &wrap<&Reader::at_time_ge>,
                arg("time"),
                return_value_policy<manage_new_object>())
        .def("at_time_lt", &wrap<&Reader::at_time_lt>,
                arg("time"),
                return_value_policy<manage_new_object>())
        .def("at_time_le", &wrap<&Reader::at_time_le>,
                arg("time"),
                return_value_policy<manage_new_object>())
        .def("grid_meta", reader_grid_meta)
        .def("grid_data", reader_grid_data)
        ;
}

