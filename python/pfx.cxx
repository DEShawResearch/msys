#include "pymod.hxx"
#include <numpy/ndarrayobject.h>
#include <pfx/pfx.hxx>

using namespace pybind11;

typedef desres::msys::pfx::Pfx pfx_t;
typedef desres::msys::pfx::Graph graph_t;

typedef struct {
    PyObject_HEAD
    pfx_t *pfx;
} PfxObject;

namespace {
    auto py_as_long = PyLong_AsLong;
    auto py_from_long = PyLong_FromLong;
}

static PyTypeObject ptype;

static PyObject* pfx_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
    PyObject *topobj, *fixobj=Py_False;
    PfxObject* self;
    graph_t* g = NULL;
    std::shared_ptr<graph_t> gptr;

    int i, n, fixbonds;

    static char* kwlist[] = {(char *)"topology", (char *)"fixbonds", 0};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|O:new", kwlist, 
                &topobj, &fixobj))
        return NULL;

    fixbonds = PyObject_IsTrue(fixobj);
    if (fixbonds<0) return NULL;

    try {
        g = cast<graph_t*>(topobj);
    } catch(cast_error) {
        if (!PySequence_Check(topobj)) {
            PyErr_Format(PyExc_TypeError, "argument must be a sequence");
            return NULL;
        }

        n = (int)PyObject_Size(topobj);
        gptr.reset(new graph_t(n));
        for (i=0; i<n; i++) {
            PyObject * sublist = PySequence_GetItem(topobj,i);
            if (!sublist || !PySequence_Check(sublist)) {
                if (sublist) {
                    Py_DECREF(sublist);
                }
                PyErr_Format(PyExc_TypeError, "input elements must be sequences");
                return NULL;
            }
            int j, nn = (int)PyObject_Size(sublist);
            for (j=0; j<nn; j++) {
                PyObject * item = PySequence_GetItem(sublist,j);
                int gid=0;
                PyErr_Clear();
                if (item) {
                    gid = (int)py_as_long(item);
                    Py_DECREF(item);
                }
                if (PyErr_Occurred()) {
                    Py_DECREF(sublist);
                    PyErr_Format(
                            PyExc_TypeError, "subsequence elements must be integer");
                    return NULL;
                }
                gptr->add_edge(i, gid); /* FIXME - error checking */
            }
            Py_DECREF(sublist);
        }
        g = gptr.get();
    }

    self = PyObject_New(PfxObject, &ptype);
    if (!self) return NULL;
    self->pfx = new pfx_t(*g, fixbonds);
    return (PyObject *)self;
}

static void pfx_dealloc(PyObject *obj) {
    delete ((PfxObject*)obj)->pfx;
    obj->ob_type->tp_free(obj);
}

PyDoc_STRVAR(apply_doc,
"apply(pos, box, vel=None) -- perform wrapping and alignment.\n"
"\n"
"box and pos are NumPy arrays with shape (3,3) and (natoms,3),\n"
"respectively, and will be modified in place.\n"
);

static PyObject* py_apply(PyObject* pySelf, PyObject* args, PyObject* kwds) {
    PyObject *boxobj=Py_None, *velobj=Py_None, *posobj=NULL;
    PyObject *boxarr=NULL, *posarr=NULL, *velarr=NULL;
    PyObject *result = NULL;
    double* box=NULL;
    pfx_t* pfx = ((PfxObject *)pySelf)->pfx;
    static char *kwlist[] = { (char *)"pos", (char *)"box", (char *)"vel", 0 };

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O|O", kwlist,
                &PyArray_Type, &posobj, &boxobj, &velobj))
        return NULL;

    /* box is either None or a 3x3 NumPy array, parsed as double */
    if (boxobj != Py_None) {
        boxarr = PyArray_FromAny(
                boxobj,
                PyArray_DescrFromType(NPY_DOUBLE), 
                2, 2, NPY_INOUT_ARRAY | NPY_FORCECAST, NULL);
        if (!boxarr) return NULL;
        if (PyArray_DIM(boxarr,0)!=3 ||
            PyArray_DIM(boxarr,1)!=3) {
            PyErr_Format(PyExc_ValueError, "box must be 3x3, got %ldx%ld",
                        PyArray_DIM(boxarr,0), PyArray_DIM(boxarr,1));
            goto error;
        }
        box = (double*)PyArray_DATA(boxarr);
    }

    /* pos is nx3 NumPyArray.  Check the type later */
    posarr = PyArray_FromAny(posobj, NULL, 2, 2, NPY_INOUT_ARRAY, NULL);
    if (!posarr) {
        Py_XDECREF(boxarr);
        return NULL;
    }
    if (PyArray_DIM(posarr,0)!=pfx->size() ||
        PyArray_DIM(posarr,1)!=3) {
        PyErr_Format(PyExc_ValueError, "pos must be %ux3, got %ldx%ld",
                pfx->size(), PyArray_DIM(posobj,0), PyArray_DIM(posobj,1));
        goto error;
    }

    /* vel is either None or an nx3 NumPy array */
    if (velobj != Py_None) {
        velarr = PyArray_FromAny(velobj, NULL, 2, 2, NPY_INOUT_ARRAY, NULL);
        if (!velarr) goto error;
        if (PyArray_DIM(velarr,0)!=pfx->size() ||
            PyArray_DIM(velarr,1)!=3) {
            PyErr_Format(PyExc_ValueError, "vel must be %ux3, got %ldx%ld",
                    pfx->size(), PyArray_DIM(velarr,0), PyArray_DIM(velarr,1));
            goto error;
        }
        if (PyArray_TYPE(velarr) != PyArray_TYPE(posarr)) {
            PyErr_Format(PyExc_ValueError, "vel must have same type as pos");
            goto error;
        }
    }

    switch (PyArray_TYPE(posarr)) {
        case NPY_FLOAT:
        {
            float* pos = (float*)PyArray_DATA(posarr);
            float* vel = velarr ? (float*)PyArray_DATA(velarr) : NULL;
            Py_BEGIN_ALLOW_THREADS
            pfx->apply(pos, box, vel);
            Py_END_ALLOW_THREADS
        }
        break;
        case NPY_DOUBLE:
        {
            double* pos = (double*)PyArray_DATA(posarr);
            double* vel = velarr ? (double*)PyArray_DATA(velarr) : NULL;
            Py_BEGIN_ALLOW_THREADS
            pfx->apply(pos, box, vel);
            Py_END_ALLOW_THREADS
        }
        break;
        default:
        {
            PyErr_Format(PyExc_ValueError, "pos must be either float or double");
            goto error;
        }
    }

    result = Py_None;
    Py_INCREF(result);

error:
    Py_XDECREF(boxarr);
    Py_XDECREF(posarr);
    Py_XDECREF(velarr);

    return result;
}

PyDoc_STRVAR(glue_doc,
"glue(atoms) -- specify atoms to be kept together during wrapping.");

static PyObject* py_glue(PyObject *self, PyObject *args, PyObject *kwds) {
    pfx_t *pfx = ((PfxObject *)self)->pfx;
    PyObject *atomobj, *atomarr;
    static char *kwlist[] = {(char *)"atoms", 0};
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &atomobj))
        return NULL;
    if (!(atomarr = PyArray_FromAny(
                    atomobj,
                    PyArray_DescrFromType(NPY_UINT32),
                    1,1,
                    NPY_C_CONTIGUOUS | NPY_ALIGNED,
                    NULL)))
        return NULL;

    int n = (int)PyArray_DIM(atomarr, 0);
    pfx->glue(n, (unsigned *)PyArray_DATA(atomarr));
    Py_DECREF(atomarr);

    Py_INCREF(Py_None);
    return Py_None;
}

PyDoc_STRVAR(align_doc,
"align(atoms, coords=None, weights=None) -- align to reference structure.\n"
"\n"
"If coords is None, then only centering will be performed.\n"
"If weights not None, it should be of size atoms.\n"
);

static PyObject* py_align(PyObject *self, PyObject *args, PyObject *kwds) {
    pfx_t *pfx = ((PfxObject *)self)->pfx;
    PyObject *atomobj, *coordobj=Py_None;
    PyObject *atomarr, *coordarr=NULL;
    PyObject *wobj=Py_None, *warr=NULL;
    const double* coords = NULL;
    const double* weights = NULL;
    static char *kwlist[] = {(char *)"atoms", (char *)"coords",
                             (char *)"weights", 0};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|OO", kwlist,
                &atomobj, &coordobj, &wobj))
        return NULL;

    if (!(atomarr = PyArray_FromAny(
                    atomobj,
                    PyArray_DescrFromType(NPY_UINT32),
                    1,1,
                    NPY_C_CONTIGUOUS | NPY_ALIGNED,
                    NULL)))
        return NULL;

    if(coordobj!=Py_None) {
        if (!(coordarr = PyArray_FromAny(
                        coordobj,
                        PyArray_DescrFromType(NPY_DOUBLE),
                        2,2,
                        NPY_C_CONTIGUOUS | NPY_ALIGNED,
                        NULL))) {
            Py_DECREF(atomarr);
            return NULL;
        }

        if (PyArray_DIM(coordarr,0)!=PyArray_DIM(atomarr,0) ||
            PyArray_DIM(coordarr,1)!=3) {
            Py_DECREF(atomarr);
            Py_DECREF(coordarr);
            PyErr_Format(PyExc_ValueError, "coords must be len(atoms) x 3");
            return NULL;
        }
        coords = (const double*)PyArray_DATA(coordarr);
    }

    if (wobj!=Py_None) {
        if (!(warr = PyArray_FromAny(
                        wobj,
                        PyArray_DescrFromType(NPY_DOUBLE),
                        1,1,
                        NPY_C_CONTIGUOUS | NPY_ALIGNED,
                        NULL))) {
            Py_DECREF(atomarr);
            Py_XDECREF(coordarr);
            return NULL;
        }
        if (PyArray_DIM(warr,0)!=PyArray_DIM(atomarr,0)) {
            Py_DECREF(atomarr);
            Py_XDECREF(coordarr);
            PyErr_Format(PyExc_ValueError, "weights must be len(atoms)");
            return NULL;
        }
        weights = (const double *)PyArray_DATA(warr);
    }

    unsigned n = (int)PyArray_DIM(atomarr, 0);
    pfx->align(n, (unsigned *)PyArray_DATA(atomarr), coords, weights);

    Py_DECREF(atomarr);
    Py_XDECREF(coordarr);
    Py_XDECREF(warr);

    Py_INCREF(Py_None);
    return Py_None;
}

PyDoc_STRVAR(rmsd_doc,
"rmsd(pos) -> root mean square distance from reference structure.\n"
"\n"
"If no reference structure has been provided with align(),\n"
"-1 is returned.\n");

static PyObject* py_rmsd(PyObject *self, PyObject *args, PyObject *kwds) {
    pfx_t *pfx = ((PfxObject *)self)->pfx;
    PyObject *posobj, *posarr;
    static char *kwlist[] = {(char *)"pos", 0};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &posobj))
        return NULL;

    if (!(posarr = PyArray_FromAny(
                   posobj, NULL,
                   2,2,
                   NPY_C_CONTIGUOUS | NPY_ALIGNED | NPY_FORCECAST,
                   NULL))) {
        return NULL;
    }

    if (PyArray_DIM(posarr,0)!=pfx->size() ||
        PyArray_DIM(posarr,1)!=3) {
        PyErr_Format(PyExc_ValueError, "pos must be len(atoms) x 3");
        Py_DECREF(posarr);
        return NULL;
    }
    double rmsd=-1;
    switch (PyArray_TYPE(posarr)) {
        case NPY_FLOAT:
            rmsd = pfx->rmsd((const float *)PyArray_DATA(posarr));
            break;
        case NPY_DOUBLE:
            rmsd = pfx->rmsd((const double*)PyArray_DATA(posarr));
            break;
        default:
            PyErr_Format(PyExc_ValueError, "pos must be float or double");
            Py_DECREF(posarr);
            return NULL;
    }
    Py_DECREF(posarr);
    return PyFloat_FromDouble(rmsd);
}


PyDoc_STRVAR(pfx_doc,
"Pfx(topology, fixbonds=False) -> periodic fix object\n"
"\n"
"topology is a list of lists, one for each atom; each sublist\n"
"contains the bonded atoms for the atom in the corresponding\n"
"list position.  If fixbonds is True, wrapping will be performed\n"
"on the bonds; this is expensive.\n"
);

static PyMethodDef pfx_methods[] = {
    { "align", 
      (PyCFunction)py_align, 
      METH_VARARGS | METH_KEYWORDS,
      align_doc },
    { "apply", 
      (PyCFunction)py_apply, 
      METH_VARARGS | METH_KEYWORDS,
      apply_doc },
    { "glue", 
      (PyCFunction)py_glue, 
      METH_VARARGS | METH_KEYWORDS,
      glue_doc },
    { "rmsd", 
      (PyCFunction)py_rmsd, 
      METH_VARARGS | METH_KEYWORDS,
      rmsd_doc },
    { NULL, NULL }
};

static PyObject* py_natoms(PfxObject* obj, void *closure) {
    return py_from_long(obj->pfx->size());
}

static PyGetSetDef pfx_getset[] = {
    { (char *)"natoms", (getter)py_natoms, NULL, (char *)"Number of atoms"},
    { NULL }
};

static PyObject* extract_3x3(PyObject* Aobj) {
    PyObject* Aarr;
    if (!(Aarr = PyArray_FromAny(
                 Aobj,
                 PyArray_DescrFromType(NPY_FLOAT64),
                 2,2,
                 NPY_C_CONTIGUOUS | NPY_ALIGNED | NPY_FORCECAST,
                 NULL))) {
        return NULL;
    }
    if (PyArray_DIM(Aarr,0)!=3 ||
        PyArray_DIM(Aarr,1)!=3) {
        PyErr_Format(PyExc_ValueError, "Expected 3x3 matrix; got %ldx%ld",
                PyArray_DIM(Aarr,0),
                PyArray_DIM(Aarr,1));
        Py_DECREF(Aarr);
        return NULL;
    }
    return Aarr;
}

PyDoc_STRVAR(inverse_doc,
"inverse_3x3(A) -> Ainv\n");
static handle wrap_inverse(object obj) {
    PyObject *Aobj, *Aarr;
    Aobj = obj.ptr();
    if (!(Aarr = extract_3x3(Aobj)))
        return NULL;
    const double* src = (const double*)PyArray_DATA(Aarr);
    double dst[9];
    desres::msys::pfx::inverse_3x3(dst, src);
    npy_intp dims[2] = {3,3};
    auto Ainv = PyArray_SimpleNew(2, dims, NPY_FLOAT64);
    memcpy(PyArray_DATA(Ainv), dst, sizeof(dst));
    return Ainv;
}

PyDoc_STRVAR(svd_doc,
"svd_3x3(A) -> U, w, V\n"
"\n"
"svd_3x3 computes the singular value decomposition of the 3x3 matrix A.\n"
"The result is always calculated and returned in double precision.\n");
static handle wrap_svd(args _args, kwargs kwds) {
    static char *kwlist[] = {(char *)"A", 0};
    PyObject* Aobj, *Aarr, *U, *W, *V, *result;
    double u[9], w[3], v[9];
    npy_intp mdims[2] = {3,3};
    npy_intp vdims[1] = {3};

    if (!PyArg_ParseTupleAndKeywords(_args.ptr(), kwds.ptr(), "O", kwlist, &Aobj))
        return NULL;
    if (!(Aarr = extract_3x3(Aobj)))
        return NULL;
    
    memcpy(u,PyArray_DATA(Aarr), sizeof(u));
    desres::msys::pfx::svd_3x3(u,w,v);
    Py_DECREF(Aarr);

    result = PyTuple_New(3);
    U = PyArray_SimpleNew(2,mdims,NPY_FLOAT64);
    W = PyArray_SimpleNew(1,vdims,NPY_FLOAT64);
    V = PyArray_SimpleNew(2,mdims,NPY_FLOAT64);
    memcpy(PyArray_DATA(U), u, sizeof(u));
    memcpy(PyArray_DATA(W), w, sizeof(w));
    memcpy(PyArray_DATA(V), v, sizeof(v));
    PyTuple_SET_ITEM(result,0,U);
    PyTuple_SET_ITEM(result,1,W);
    PyTuple_SET_ITEM(result,2,V);
    return result;
}


PyDoc_STRVAR(aligned_rmsd_doc,
"aligned_rmsd(X, Y, weight=None) -> mat, rmsd\n\n"
"Compute the matrix aligning Y onto X, optionally with weights.  Return the matrix and the rmsd.\n"
);

static 
handle wrap_aligned_rmsd(args _args, kwargs kwds) {
    static char *kwlist[] = {(char *)"X", (char *)"Y", (char *)"weight", 0};
    PyObject *Xobj, *Yobj, *Wobj=NULL;
    PyObject *Xarr, *Yarr, *Marr, *Warr=NULL;
    if (!PyArg_ParseTupleAndKeywords(_args.ptr(), kwds.ptr(), "OO|O", kwlist, &Xobj, &Yobj, &Wobj))
        return NULL;
    if (!(Xarr = PyArray_FromAny(
                    Xobj,
                    PyArray_DescrFromType(NPY_FLOAT64),
                    2,2,
                    NPY_C_CONTIGUOUS | NPY_ALIGNED | NPY_FORCECAST,
                    NULL)))
        return NULL;
    if (!(Yarr = PyArray_FromAny(
                    Yobj,
                    PyArray_DescrFromType(NPY_FLOAT64),
                    2,2,
                    NPY_C_CONTIGUOUS | NPY_ALIGNED | NPY_FORCECAST,
                    NULL))) {
        Py_DECREF(Xarr);
        return NULL;
    }
    if (Wobj && !(Warr = PyArray_FromAny(
                    Wobj,
                    PyArray_DescrFromType(NPY_FLOAT64),
                    1,1,
                    NPY_C_CONTIGUOUS | NPY_ALIGNED | NPY_FORCECAST,
                    NULL))) {
        Py_DECREF(Xarr);
        Py_DECREF(Yarr);
        return NULL;
    }
    unsigned n = PyArray_DIM(Xarr,0);
    if (PyArray_DIM(Xarr,1)!=3 || PyArray_DIM(Yarr,1)!=3 ||
        n==0 || PyArray_DIM(Yarr,0)!=n) {
        PyErr_Format(PyExc_ValueError, "Require equal sized nx3 matrices");
        Py_DECREF(Xarr);
        Py_DECREF(Yarr);
        Py_XDECREF(Warr);
        return NULL;
    }
    if (Warr && (PyArray_DIM(Warr,0)!=n)) {
        PyErr_Format(PyExc_ValueError, "Require weights of same length as X, Y");
        Py_DECREF(Xarr);
        Py_DECREF(Yarr);
        Py_XDECREF(Warr);
        return NULL;
    }

    npy_intp dims[2] = {3,3};
    Marr = PyArray_SimpleNew(2,dims,NPY_FLOAT64);
    const double* X = (double *)PyArray_DATA(Xarr);
    const double* Y = (double *)PyArray_DATA(Yarr);
    const double* W = Warr ? (double *)PyArray_DATA(Warr) : (double *)NULL;

    double* mat = (double *)PyArray_DATA(Marr);
    double rmsd = desres::msys::pfx::compute_alignment(n, NULL, X, Y, mat, W);
    Py_DECREF(Xarr);
    Py_DECREF(Yarr);
    Py_XDECREF(Warr);
    PyObject* obj = Py_BuildValue("(O,d)", Marr, rmsd);
    Py_DECREF(Marr);
    return obj;
}

PyDoc_STRVAR(module_doc,
"A high level interface for wrapping, centering, and alignment.\n"
"\n"
"This module provides simple, yet performant methods for manipulating\n"
"trajectories of systems with connected atoms and periodic boundary\n"
"conditions.");


PYBIND11_MODULE(pfx, m) {
    _import_array();
    if (PyErr_Occurred()) throw error_already_set();

    Py_TYPE(&ptype) = &PyType_Type;
    ptype.tp_name = "pfx.Pfx";
    ptype.tp_doc = pfx_doc;
    ptype.tp_basicsize = sizeof(PfxObject);
    ptype.tp_alloc = PyType_GenericAlloc;
    ptype.tp_new = pfx_new;
    ptype.tp_dealloc = pfx_dealloc;
    ptype.tp_flags = Py_TPFLAGS_DEFAULT;
    ptype.tp_getset = pfx_getset;
    ptype.tp_methods = pfx_methods;
    if (PyType_Ready(&ptype)) throw error_already_set();

    Py_INCREF((PyObject *)&ptype);
    PyModule_AddObject(m.ptr(), "Pfx", (PyObject *)&ptype);
    m.def("svd_3x3", wrap_svd, svd_doc);
    m.def("inverse_3x3", wrap_inverse, inverse_doc);
    m.def("aligned_rmsd", wrap_aligned_rmsd, aligned_rmsd_doc);
    m.attr("__doc__") = module_doc;
}

