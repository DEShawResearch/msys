#define NO_IMPORT_ARRAY 1
#include "unique_symbol.hxx"
#include "spatial_hash.hxx"
#include "wrap_obj.hxx"

using namespace boost::python;
using namespace desres::msys;

typedef desres::msys::SpatialHash<float> Hash;

static Hash* hash_new(PyObject* posobj, PyObject* idobj) {

    PyObject* posarr = PyArray_FromAny(posobj,
                PyArray_DescrFromType(NPY_FLOAT32),
                2, 2, NPY_C_CONTIGUOUS | NPY_ALIGNED,
                NULL);
    if (!posarr) throw_error_already_set();

    PyObject* idarr = PyArray_FromAny(idobj,
                PyArray_DescrFromType(NPY_UINT32),
                1, 1, NPY_C_CONTIGUOUS | NPY_ALIGNED,
                NULL);
    if (!idarr) {
        Py_DECREF(posarr);
        throw_error_already_set();
    }

    const float* pos = (const float*)PyArray_DATA(posarr);
    const Id* ids = (const Id*)PyArray_DATA(idarr);
    Id n = PyArray_DIM(idarr,0);

    Hash* hash = new Hash(pos, n, ids);

    Py_DECREF(idarr);
    Py_DECREF(posarr);
    return hash;
}


template <typename Func, typename Param>
static PyObject* hash_find(Hash& hash,
                           Param r,
                           PyObject* posobj,
                           PyObject* idsobj,
                           PyObject* boxobj,
                           Func func) {

    PyObject* posarr = PyArray_FromAny(posobj,
                PyArray_DescrFromType(NPY_FLOAT32),
                2, 2, NPY_C_CONTIGUOUS | NPY_ALIGNED,
                NULL);
    if (!posarr) throw_error_already_set();

    PyObject* idsarr = PyArray_FromAny(idsobj,
                PyArray_DescrFromType(NPY_UINT32),
                1, 1, NPY_C_CONTIGUOUS | NPY_ALIGNED,
                NULL);
    if (!idsarr) {
        Py_DECREF(posarr);
        throw_error_already_set();
    }

    PyObject* boxarr = NULL;
    if (boxobj!=Py_None) {
        boxarr = PyArray_FromAny(boxobj,
                 PyArray_DescrFromType(NPY_FLOAT64),
                 2, 2, NPY_C_CONTIGUOUS | NPY_ALIGNED,
                 NULL);
        if (!boxarr) {
            Py_DECREF(idsarr);
            Py_DECREF(posarr);
            throw_error_already_set();
        }
    }

    const float* pos = (const float*)PyArray_DATA(posarr);
    const double* box = boxarr ? (const double *)PyArray_DATA(boxarr) : NULL;
    const Id* ids = (const Id*)PyArray_DATA(idsarr);
    Id n = PyArray_DIM(idsarr,0);

    IdList result = (hash.*func)(r, pos, n, ids, box);

    Py_XDECREF(boxarr);
    Py_DECREF(idsarr);
    Py_DECREF(posarr);
    npy_intp dims[1];
    dims[0] = result.size();
    PyObject* arr = PyArray_SimpleNew(1,dims,NPY_UINT32);
    if (!arr) throw_error_already_set();
    if (!result.empty()) memcpy(PyArray_DATA(arr), 
                                &result[0],
                                result.size()*sizeof(Id));
    return arr;
}

static PyObject* hash_find_within(Hash& hash,
                                  float r,
                                  PyObject* pos,
                                  PyObject* ids,
                                  PyObject* box,
                                  bool reuse_voxels) {

    return reuse_voxels ? hash_find(hash, r, pos, ids, box, &Hash::find_within)
                        : hash_find(hash, r, pos, ids, box, &Hash::findWithin);
}

static PyObject* hash_find_nearest(Hash& hash,
                                   int k,
                                   PyObject* pos,
                                   PyObject* ids,
                                   PyObject* box) {

    return hash_find(hash, k, pos, ids, box, &Hash::findNearest);
}

namespace desres { namespace msys {
    void export_spatial_hash() {

        class_<Hash>("SpatialHash", no_init)
            .def("__init__", make_constructor(
                        hash_new, default_call_policies(),
                        (arg("pos"), arg("ids"))))
            .def("voxelize", &Hash::voxelize, return_internal_reference<>())
            .def("findWithin", hash_find_within,
                    (arg("r"), 
                     arg("pos"), 
                     arg("ids"), 
                     arg("box")=object(),
                     arg("reuse_voxels")=false))
            .def("findNearest", hash_find_nearest,
                    (arg("k"), 
                     arg("pos"), 
                     arg("ids"), 
                     arg("box")=object()))
            ;
    }
}}

