#include <boost/python.hpp>
#include <numpy/ndarrayobject.h>
#include "cealign.hxx"

using namespace desres::msys;
using namespace boost::python;

typedef CEAlign<double> Engine;
typedef std::shared_ptr<PyObject> objptr;
static void dtor(PyObject* obj) { Py_DECREF(obj); }

static list compute(Engine const& engine, list Aatoms, object Aobj,
                                          list Batoms, object Bobj) {
    objptr Aarr, Barr;
    Aarr.reset(PyArray_FromAny(
                Aobj.ptr(),
                PyArray_DescrFromType(NPY_FLOAT64),
                2,2,NPY_C_CONTIGUOUS, NULL),
            dtor);
    if (!Aarr) throw_error_already_set();
    Barr.reset(PyArray_FromAny(
                Bobj.ptr(),
                PyArray_DescrFromType(NPY_FLOAT64),
                2,2,NPY_C_CONTIGUOUS, NULL),
            dtor);
    if (!Barr) throw_error_already_set();
    if (PyArray_DIM(Aarr.get(), 1)!=3 ||
        PyArray_DIM(Barr.get(), 1)!=3) {
        PyErr_Format(PyExc_ValueError, "Expected nx3 positions");
        throw_error_already_set();
    }

    IdList A(len(Aatoms)), B(len(Batoms));
    for (boost::python::ssize_t i=0, n=len(Aatoms); i<n; i++) A[i] = extract<Id>(Aatoms[i]);
    for (boost::python::ssize_t i=0, n=len(Batoms); i<n; i++) B[i] = extract<Id>(Batoms[i]);
    if (*std::max_element(A.begin(), A.end()) >= PyArray_DIM(Aarr.get(), 0)) {
        PyErr_Format(PyExc_ValueError, "Out of range index in A");
        throw_error_already_set();
    }
    if (*std::max_element(B.begin(), B.end()) >= PyArray_DIM(Barr.get(), 0)) {
        PyErr_Format(PyExc_ValueError, "Out of range index in B");
        throw_error_already_set();
    }
    if (A.size()<engine.window_size() || B.size()<engine.window_size()) {
        return list();
    }

    const double* Apos = (double *)PyArray_DATA(Aarr.get());
    const double* Bpos = (double *)PyArray_DATA(Barr.get());
    Engine::MapList maps = engine.compute(A, Apos, B, Bpos);

    list result;
    for (size_t i=0, n=maps.size(); i<n; i++) {
        Engine::Map const& map = maps[i];
        list aids, bids, elem;
        for (size_t j=0, m=map.first.size(); j<m; j++) {
            aids.append(map.first[j]);
            bids.append(map.second[j]);
        }
        elem.append(aids);
        elem.append(bids);
        result.append(elem);
    }
    return result;
}


BOOST_PYTHON_MODULE(_cealign) {

    import_array();
    if (PyErr_Occurred()) return;

    class_<Engine>("CEAlign", no_init)
        .def("WithDefaults", &Engine::WithDefaults).staticmethod("WithDefaults")
        .def("compute", compute)
        ;
}

