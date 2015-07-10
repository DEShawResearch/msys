#include "molfilemodule.hxx"
#include "numpy/arrayobject.h"
#include <map>
#include <string>
#include <sstream>
#include <dlfcn.h>

#include "molfile/findframe.hxx"
#include "molfile/dtrplugin.hxx"

#include <boost/python.hpp>
#include <boost/thread.hpp>

using namespace desres::molfile;
using namespace boost::python;
namespace ff=findframe;

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

    DtrWriter* dtr_init(std::string const& path,
                        uint32_t natoms,
                        int mode,
                        uint32_t fpf) {
        DtrWriter* w = new DtrWriter(natoms, fpf);
        try {
            w->init(path, (DtrWriter::Mode)mode);
        }
        catch (std::exception& e) {
            PyErr_Format(PyExc_IOError, "Unable to initialize DtrWriter at %s: %s", path.c_str(), e.what());
            delete w;
            throw error_already_set();
        }
        return w;
    }

    void dtr_append(DtrWriter& w, double time, dict keyvals) {
        dtr::KeyMap map;
        list items = keyvals.items();
        for (unsigned i=0, n=len(items); i<n; i++) {
            object item = items[i];
            std::string key = extract<std::string>(item[0]);
            object val = item[1];
            PyObject* ptr = val.ptr();
            dtr::Key keyval;
            if (PyString_Check(ptr)) {
                const char* s = PyString_AsString(ptr);
                keyval.set(PyString_AsString(ptr), strlen(s));
            } else if (PyArray_Check(ptr)) {
                if (PyArray_NDIM(ptr)!=1) {
                    PyErr_Format(PyExc_ValueError, "array for key %s must be 1d", key.c_str());
                    throw_error_already_set();
                }
                int npytype = PyArray_TYPE(ptr);
                const void* data = PyArray_DATA(ptr);
                unsigned n = PyArray_DIM(ptr,0);
                switch (npytype) {
                    case NPY_INT32:
                        keyval.set((const int32_t*)data, n);
                        break;
                    case NPY_UINT32:
                        keyval.set((const uint32_t*)data, n);
                        break;
                    case NPY_INT64:
                        keyval.set((const int64_t*)data, n);
                        break;
                    case NPY_UINT64:
                        keyval.set((const uint64_t*)data, n);
                        break;
                    case NPY_FLOAT32:
                        keyval.set((const float*)data, n);
                        break;
                    case NPY_FLOAT64:
                        keyval.set((const double*)data, n);
                        break;
                    default:
                        PyErr_Format(PyExc_ValueError, "unsupported array type for key %s", key.c_str());
                        throw_error_already_set();
                }
            } else {
                PyErr_Format(PyExc_ValueError, "values must be string or numpy array");
                throw error_already_set();
            }
            map[key] = keyval;
        }
        w.append(time, map);
    }

    void export_dtrwriter() {

        class_<DtrWriter, boost::noncopyable>("DtrWriter", no_init)
            .def("__init__",
                make_constructor( 
                    dtr_init,
                    default_call_policies(),
                    (arg("path"),
                     arg("natoms"),
                     arg("mode")=0,
                     arg("frames_per_file")=object())))
            .def("append", dtr_append,
                    (arg("time"),
                     arg("keyvals")))
            .def("sync", &DtrWriter::sync)
            ;

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
    export_dtrwriter();

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

namespace {

    FrameSetReader * dtr_from_path( object& pathobj, bool sequential ) {
        std::string path = extract<std::string>(pathobj);
        FrameSetReader * reader = NULL;
        if (StkReader::recognizes(path)) {
            reader = new StkReader(path,
                                sequential ? DtrReader::SequentialAccess : 0);
        } else {
            reader = new DtrReader(path, 
                                sequential ? DtrReader::SequentialAccess : 0);
        }
        try {
            reader->init();
        }
        catch (std::exception &e) {
            delete reader;
            throw;
        }
        return reader;
    }

    FrameSetReader * dtrreader_init( object& path, bool sequential ) {
        if (!path.is_none()) {
	    return dtr_from_path(path, sequential);
	} else {
            PyErr_Format(PyExc_ValueError, "Must supply path");
            throw error_already_set();
        }
        return NULL;
    }

    const char * fileinfo_doc = 
        "fileinfo(index) -> path, time, offset, framesize, first, last, filesize, dtrpath, dtrsize\n"
        "  file contains frames [first, last)\n";

    tuple fileinfo(FrameSetReader& self, Py_ssize_t index) {
        const DtrReader *comp = self.component(index);
        if (!comp) {
            PyErr_SetString(PyExc_IndexError, "index out of bounds");
            throw error_already_set();
        }
        const key_record_t &key = comp->keys[index];
        std::string path = comp->framefile(index);

        /* remainder: index of this frame in the framefile */
        Py_ssize_t remainder = index % comp->framesperfile();

        /* first: the index of the first frame in the file */
        Py_ssize_t first = index - remainder;

        /* last: the index of the last frame in the file */
        Py_ssize_t last = first;

        /* in the DtrReader class, the timekeys array is truncated if the
         * last few times overlap later dtrs.  However, we don't want overlap
         * to affect the reported size of the file in which this frame
         * residues.  We therefore check to see if overlapping timekeys
         * need to be considered in computing the size of the file. */

        Py_ssize_t dtrsize = comp->keys.full_size();

        Py_ssize_t filesize;
        for (last=index+1; last<dtrsize; ++last) {
            if (comp->framefile(last) != path) break;
        }
        --last;
        if (last==index) {
            /* file size is that of a single frame */
            filesize = key.offset() + key.size();

        } else {
            /* filesize from key in frames */
            filesize = comp->keys[last].offset() + comp->keys[last].size();

        }

        const std::string &dtrpath = comp->path();

        return make_tuple(
                path.c_str(), key.time(), key.offset(), key.size(),
                first, last+1, filesize, dtrpath.c_str(), dtrsize );
    }

    const char * frame_doc = 
        "frame(index, bytes=None, with_gids=False, keyvals=None) -> Frame\n"
        "Read bytes from disk if bytes are not provided\n"
        "If keyvals is not None, it should be a dict, and raw data from\n"
        "the frame will be provided.\n";

    object frame(FrameSetReader& self, Py_ssize_t index, object& bytes,
            bool with_gids, object keyvals) {
        /* The component method modifies index.  What a dumb API. */
        Py_ssize_t global_index = index;
        const DtrReader *comp = self.component(index);
        if (!comp) {
            PyErr_SetString(PyExc_IndexError, "index out of bounds");
            throw error_already_set();
        }
        Frame *frame = new Frame(comp->natoms(), self.has_velocities(),
                                 false, /* double precision */
                                 with_gids);
        dtr::KeyMap keymap;
        void* keybuf = NULL;
        void** keybufptr = keyvals.ptr() == Py_None ? NULL : &keybuf;
        boost::shared_ptr<void> keybuf_dtor(keybuf, free);
        if (bytes.is_none()) {
            PyThreadState *_save;
            _save = PyEval_SaveThread();
            try {
                keymap = comp->frame(index, *frame, keybufptr);
            }
            catch (std::exception &e) {
                PyEval_RestoreThread(_save);
                delete frame;
                PyErr_Format(PyExc_IOError, "Error reading frame: global index %ld dtr path %s local index %ld frame file %s\n%s",
                        global_index, comp->path().c_str(),
                        index, comp->framefile(index).c_str(),
                        e.what());
                throw_error_already_set();
            }
            PyEval_RestoreThread(_save);

        } else {
            Py_ssize_t size;
            char * data;
            if (PyString_AsStringAndSize(bytes.ptr(), &data, &size)) {
                delete frame;
                throw error_already_set();
            }
            try {
                comp->frame_from_bytes(data, size, *frame);
                frame->setTime(comp->keys[index].time());
            }
            catch (std::exception &e) {
                delete frame;
                PyErr_Format(PyExc_RuntimeError, "Failed parsing frame data of size %ld: %s", size, e.what());
                throw_error_already_set();
            }
        }
        for (dtr::KeyMap::const_iterator it=keymap.begin(), e=keymap.end(); it!=e; ++it) {
            dtr::Key const& val = it->second;
            npy_intp dims = val.count;
            PyObject* arr = NULL;
            switch (val.type) {
                case dtr::Key::TYPE_INT32:
                    arr = PyArray_SimpleNew(1, &dims, NPY_INT32);
                    val.get((int32_t*)PyArray_DATA(arr)); 
                    break;
                case dtr::Key::TYPE_UINT32:
                    arr = PyArray_SimpleNew(1, &dims, NPY_UINT32);
                    val.get((uint32_t*)PyArray_DATA(arr)); 
                    break;
                case dtr::Key::TYPE_INT64:
                    arr = PyArray_SimpleNew(1, &dims, NPY_INT64);
                    val.get((int64_t*)PyArray_DATA(arr)); 
                    break;
                case dtr::Key::TYPE_UINT64:
                    arr = PyArray_SimpleNew(1, &dims, NPY_UINT64);
                    val.get((uint64_t*)PyArray_DATA(arr)); 
                    break;
                case dtr::Key::TYPE_FLOAT32:
                    arr = PyArray_SimpleNew(1, &dims, NPY_FLOAT32);
                    val.get((float*)PyArray_DATA(arr)); 
                    break;
                case dtr::Key::TYPE_FLOAT64:
                    arr = PyArray_SimpleNew(1, &dims, NPY_FLOAT64);
                    val.get((double*)PyArray_DATA(arr)); 
                    break;
                case dtr::Key::TYPE_CHAR:
                case dtr::Key::TYPE_UCHAR:
                    arr = PyString_FromString(val.toString().c_str());
                    break;
                default:;
            }
            if (!arr) continue;
            char* key = const_cast<char*>(it->first.c_str());
            int rc = PyMapping_SetItemString(keyvals.ptr(), key, arr);
            Py_DECREF(arr);
            if (rc==-1) throw_error_already_set();
        }
        return object(frame);
    }

    const char reload_doc[] =
        "reload() -> number of timekeys reloaded -- reload frames in the dtr/stk";
    int reload(FrameSetReader& self) {
        int changed=0;
        self.init(&changed);
        return changed;
    }

    object get_times(const FrameSetReader& self) {
        Py_ssize_t n = self.size();
        Py_ssize_t dims[1] = { n };
        PyObject * arr = backed_vector( 1, dims, DOUBLE, NULL, NULL );
        if (self.times(0,n,(double *)array_data(arr))!=n) {
            Py_DECREF(arr);
            arr=NULL;
            PyErr_Format(PyExc_RuntimeError, "Error reading times");
            throw error_already_set();
        }
        return object(handle<>(arr));
    }
    Py_ssize_t frameset_size(const FrameSetReader& self, Py_ssize_t n) {
        return self.frameset(n)->size();
    }
    std::string frameset_path(const FrameSetReader& self, Py_ssize_t n) {
        return self.frameset(n)->path();
    }
    bool frameset_is_compact(const FrameSetReader& self, Py_ssize_t n) {
        return self.frameset(n)->keys.is_compact();
    }

    struct Oracle {
        const FrameSetReader& self;
        Oracle( const FrameSetReader& f ) : self(f) {}
        double operator[](ssize_t i) const {
            double x;
            self.times(i,1,&x);
            return x;
        }
    };

    ssize_t index_near( const FrameSetReader& self, double T ) {
        return ff::at_time_near(self.size(), Oracle(self), T);
    }
    ssize_t index_le( const FrameSetReader& self, double T ) {
        return ff::at_time_le(self.size(), Oracle(self), T);
    }
    ssize_t index_lt( const FrameSetReader& self, double T ) {
        return ff::at_time_lt(self.size(), Oracle(self), T);
    }
    ssize_t index_ge( const FrameSetReader& self, double T ) {
        return ff::at_time_ge(self.size(), Oracle(self), T);
    }
    ssize_t index_gt( const FrameSetReader& self, double T ) {
        return ff::at_time_gt(self.size(), Oracle(self), T);
    }
    ssize_t total_bytes( const FrameSetReader& self ) {
	return (self.total_bytes());
    }

    std::string my_path(FrameSetReader const& self) {
        return self.path();
    }

    FrameSetReader* from_timekeys(std::string const& path,
                                  list lfnames, list ltimekeys) {
        if (len(lfnames) != len(ltimekeys)) {
            PyErr_Format(PyExc_ValueError, "fnames and timekeys must be the same length");
            throw_error_already_set();
        }
        unsigned i,n = len(lfnames);
        std::vector<std::string> fnames(n);
        std::vector<Timekeys> timekeys(n);
        for (i=0; i<n; i++) {
            fnames[i] = extract<std::string>(lfnames[i]);
            timekeys[i] = extract<Timekeys>(ltimekeys[i]);
        }
        StkReader* stk = new StkReader(path);
        stk->append(fnames, timekeys);
        return stk;
    }

    void tk_init(Timekeys& tk, std::string const& path) {
        PyThreadState *_save;
        _save = PyEval_SaveThread();
        boost::shared_ptr<PyThreadState> rst(_save, PyEval_RestoreThread);
        tk.init(path);
    }

}

void desres::molfile::export_dtrreader() {

    class_<Timekeys>("Timekeys", init<>())
        .def("init", tk_init)
        .def("size",                    &Timekeys::size)
        .add_property("framesperfile",  &Timekeys::framesperfile)
        .add_property("framesize",      &Timekeys::framesize)
        .add_property("interval",       &Timekeys::interval)
        ;

    class_<FrameSetReader, boost::noncopyable>("DtrReader", no_init)
        .def("__init__", 
                make_constructor( 
                    dtrreader_init,
                    default_call_policies(),
                    (arg("path")=object(), 
                     /* WARNING: use sequential=True only from one thread! */
                     arg("sequential")=false )))
        .add_property("path", my_path)
        .add_property("natoms", &FrameSetReader::natoms)
        .add_property("nframes", &FrameSetReader::size)
        .add_property("nframesets", &FrameSetReader::nframesets)

        .def("fromTimekeys", from_timekeys, 
                return_value_policy<manage_new_object>())
        .staticmethod("fromTimekeys")

        .def("total_bytes", total_bytes)
        .def("index_near", index_near)
        .def("index_le", index_le)
        .def("index_lt", index_lt)
        .def("index_ge", index_ge)
        .def("index_gt", index_gt)
        .def("frameset_size", frameset_size)
        .def("frameset_path", frameset_path)
        .def("frameset_is_compact", frameset_is_compact)
        .def("fileinfo", fileinfo, fileinfo_doc)
        .def("frame", frame, frame_doc,
                (arg("index") 
                ,arg("bytes")=object()
                ,arg("with_gids")=false
                ,arg("keyvals")=object()))
        .def("reload", reload, reload_doc)
        .def("times", get_times)
        ;

                
    scope().attr("dtr_serialized_version")=dtr_serialized_version();
}

