#include "molfilemodule.hxx"
#include "numpy/arrayobject.h"
#include <map>
#include <string>
#include <dlfcn.h>

#include <boost/python.hpp>
#include <boost/thread.hpp>

using namespace desres::molfile;
using namespace boost::python;

using boost::python::ssize_t;

// all the numpy stuff has to live in this file!

PyObject **desres::molfile::object_array(int size, PyObject **returned_result) {
    npy_intp dims=size;
    PyObject *result = PyArray_SimpleNew(1, &dims, NPY_OBJECT);
    if (!result) return NULL;
    *returned_result = result;
    return reinterpret_cast<PyObject**>(PyArray_DATA(result));
}

PyObject *desres::molfile::backed_vector( 
        int nd, Py_ssize_t *dims, DataType type, void *data, PyObject *base ) {
    npy_intp mydims[NPY_MAXDIMS];
    for (int i=0; i<nd; i++) mydims[i] = dims[i];
    int mytype = type==INT ? PyArray_INT :
        type==FLOAT ? PyArray_FLOAT :
        type==DOUBLE ? PyArray_DOUBLE :
        0;
    if (!mytype) {
        PyErr_Format(PyExc_RuntimeError, "Unsupported type '%d'", type);
        return NULL;
    }
    PyObject *result;
    if (base) { 
        result = PyArray_SimpleNewFromData( nd, mydims, mytype, (char *)data);
        Py_INCREF( PyArray_BASE(result) = reinterpret_cast<PyObject*>(base) );
    } else {
        result = PyArray_SimpleNewFromData( nd, mydims, mytype, NULL);
        if (data) {
            memcpy( PyArray_DATA(result), data, PyArray_NBYTES(result) );
        }
    }
    return result;
}

void * desres::molfile::array_data( PyObject * arr ) {
    return PyArray_DATA(arr);
}

namespace {

    // create a Plugin object representing each molfile_plugin_t.
    int register_cb( void *v, vmdplugin_t *raw_plugin ) {
        object& func = *(object *)v;
        func(reinterpret_cast<molfile_plugin_t*>(raw_plugin));
        return MOLFILE_SUCCESS;
    }

    typedef int (*initfunc)(void);
    typedef int (*regfunc)(void *, vmdplugin_register_cb);
    typedef int (*finifunc)(void);

    // scan the .so file at the given path and load 
    void register_shared_library(const std::string& path, object& callback) {
        void * handle = dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
        if (!handle) throw std::runtime_error(dlerror());
        void * ifunc = dlsym(handle, "vmdplugin_init");
        if (ifunc && ((initfunc)(ifunc))()) {
            dlclose(handle);
            throw std::runtime_error("vmdplugin_init() failed");
        }
        void * rfunc = dlsym(handle, "vmdplugin_register");
        if (!rfunc) {
            dlclose(handle);
            throw std::runtime_error("No vmdplugin_register() found");
        }
        ((regfunc)rfunc)(&callback, register_cb);
    }

    /* return a list of all the plugins */
    void register_all(object& func) {
        MOLFILE_REGISTER_ALL(&func, register_cb); 
    }

    void destructor(PyObject* obj) { Py_XDECREF(obj); }
    typedef boost::shared_ptr<PyObject> objptr;

#define DIVVY_BEGIN_END(ntasks,rank,nprocs,begin,end) do {\
  unsigned int __quo__, __rem__;                                          \
  unsigned int __nt__=(ntasks),__rk__=(rank),__np__=(nprocs);             \
  __quo__=__nt__/__np__, __rem__=__nt__%__np__;                           \
  (begin)= __quo__*__rk__ + (__rk__<__rem__?__rk__:__rem__);              \
  (end)  = __quo__*(__rk__+1) + (__rk__<__rem__?(__rk__+1):__rem__);      \
} while(0)

    // read a single frame from the given reader and store in provided buffers
    static void read_frame(Reader& r, int64_t fid, 
            unsigned ngids, const unsigned* gids, float* pos, double* box) {
        boost::shared_ptr<Frame> frame(r.frame(fid));
        memcpy(box, frame->box(), 9*sizeof(double));
        const float* src = frame->pos();
        for (unsigned i=0; i<ngids; i++) {
            unsigned gid = gids[i];
            if (gid>=frame->natoms()) throw std::runtime_error("invalid gid");
            memcpy(pos, src+3*gid, 3*sizeof(float));
            pos += 3;
        }
    }
    
    // Worker myrank reads a contiguous range of frames within array fids
    // and stores in corresponding elements of posptrs and boxptrs
    void worker(unsigned nworkers, unsigned myrank, size_t nfids, Reader* r, 
                const int64_t* fids, unsigned ngids, const unsigned* gids,
                const std::vector<float*>*  posptrs,
                const std::vector<double*>* boxptrs) {
        unsigned begin, end;
        DIVVY_BEGIN_END(nfids, myrank, nworkers, begin, end);
        for (; begin!=end; ++begin) {
            read_frame(*r, fids[begin], ngids, gids, 
                    (*posptrs)[begin], (*boxptrs)[begin]);
        }
    }

    // Read frames from reader, storing in provided buffers.
    // fidobj: list of frame ids, nfids long
    // gidobj: list of atom ids, ngids long
    // posbuffers: list of nfids float buffers of size ngidsx3
    // boxbuffers: list of nfids double buffers of size 3x3
    // maxthreads: 0 to read in current thread, otherwise use 
    //             min(maxthreads, nfids to read frames)
    void read_frames(Reader& r, PyObject* fidobj, PyObject* gidobj,
                     PyObject* posbuffers,
                     PyObject* boxbuffers,
                     ssize_t maxthreads) {

        objptr fidarr, gidarr;
        fidarr.reset(PyArray_FromAny(fidobj, PyArray_DescrFromType(NPY_INT64),
                                     1, 1, NPY_C_CONTIGUOUS, NULL), destructor);
        if (!fidarr) throw_error_already_set();
        const int64_t* fids = (const int64_t*)PyArray_DATA(fidarr.get());
        gidarr.reset(PyArray_FromAny(gidobj, PyArray_DescrFromType(NPY_UINT32),
                                     1, 1, NPY_C_CONTIGUOUS | NPY_FORCECAST, 
                                     NULL), destructor);
        if (!gidarr) throw_error_already_set();
        const unsigned* gids = (const unsigned*)PyArray_DATA(gidarr.get());
        ssize_t nfids = PyArray_DIM(fidarr.get(), 0);
        unsigned ngids = PyArray_DIM(gidarr.get(), 0);
        std::vector<float*> posptrs(nfids);
        std::vector<double*> boxptrs(nfids);
        for (ssize_t i=0; i<nfids; i++) {
            objptr buf(PySequence_GetItem(posbuffers, i), destructor);
            if (!buf) throw_error_already_set();
            if (!PyArray_Check(buf.get()) ||
                 PyArray_TYPE(buf.get())!=NPY_FLOAT ||
                 PyArray_NDIM(buf.get())!=2 ||
                 PyArray_DIM(buf.get(),0)!=ngids ||
                 PyArray_DIM(buf.get(),1)!=3
                 ) {
                PyErr_Format(PyExc_ValueError, "Expected %ux3 numpy float array for position index %ld", ngids, i);
                throw_error_already_set();
            }
            posptrs[i] = (float *)PyArray_DATA(buf.get());
        }
        for (ssize_t i=0; i<nfids; i++) {
            objptr buf(PySequence_GetItem(boxbuffers, i), destructor);
            if (!buf) throw_error_already_set();
            if (!PyArray_Check(buf.get()) ||
                 PyArray_TYPE(buf.get())!=NPY_DOUBLE ||
                 PyArray_NDIM(buf.get())!=2 ||
                 PyArray_DIM(buf.get(),0)!=3 ||
                 PyArray_DIM(buf.get(),1)!=3
                 ) {
                PyErr_Format(PyExc_ValueError, "Expected 3x3 numpy double array for box index %ld", i);
                throw_error_already_set();
            }
            boxptrs[i] = (double *)PyArray_DATA(buf.get());
        }

        ssize_t nthreads = std::min(nfids, maxthreads);
        if (nthreads>0) {
            typedef boost::shared_ptr<boost::thread> ThreadPtr;
            std::vector<ThreadPtr> threads(nthreads);
            for (unsigned i=0; i<nthreads; i++) {
                threads[i].reset(new boost::thread(worker,
                            nthreads, i, nfids, 
                            &r, fids, ngids, gids, &posptrs, &boxptrs));
            }
            for (unsigned i=0; i<nthreads; i++) {
                threads[i]->join();
            }
        } else {
            worker(1, 0, nfids, &r, fids, ngids, gids, &posptrs, &boxptrs);
        }
    }

}

BOOST_PYTHON_MODULE(_molfile) {

    import_array();
    if (PyErr_Occurred()) return;

    export_plugin();
    export_reader();
    export_frame();
    export_writer();
    export_dtrreader();

    // I could use a static to initialize these types on first access,
    // but then the type objects wouldn't be present in the module
    // until first use.  This would cause help(molfile) to be incomplete.
    if (initialize_atom()) return;

    static PyObject *pyAtom = reinterpret_cast<PyObject *>(&AtomType);
    Py_INCREF(pyAtom);
    scope().attr("Atom")=object(handle<>(pyAtom));

    MOLFILE_INIT_ALL;
    def("register_all", register_all);
    def("register_shared_library", register_shared_library);
    def("read_frames", read_frames);
}

