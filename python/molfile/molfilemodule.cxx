#if _MSC_VER
#define BOOST_NO_CXX11_TEMPLATE_ALIASES
#endif

#include "molfilemodule.hxx"
#include "numpy/arrayobject.h"
#include <map>
#include <string>
#include <sstream>
#include <boost/python.hpp>

#ifdef WIN32
#define ssize_t boost::python::ssize_t
#endif

#include "molfile/findframe.hxx"
#include "molfile/dtrplugin.hxx"
#include "molfile/dtrframe.hxx"

#ifdef WIN32
#undef ssize_t
#endif
#include <boost/thread.hpp>
#include <boost/python/docstring_options.hpp>

using namespace desres::molfile;
using namespace boost::python;
namespace ff=findframe;

using boost::python::ssize_t;

namespace {
#if PY_MAJOR_VERSION >= 3
    auto py_string_check = [](PyObject* o) { return PyUnicode_Check(o); };
    auto py_bytes_from_string_size = PyBytes_FromStringAndSize;
    auto py_string_from_string_size = PyUnicode_FromStringAndSize;
    auto py_as_string = PyUnicode_AsUTF8;
    auto py_as_string_size = PyBytes_AsStringAndSize;
    auto py_as_bytes = PyBytes_FromStringAndSize;
#else
    auto py_string_check = [](PyObject* o) { return PyString_Check(o); };
    auto py_bytes_from_string_size = PyString_FromStringAndSize;
    auto py_string_from_string_size = PyString_FromStringAndSize;
    auto py_as_string = PyString_AsString;
    auto py_as_string_size = PyString_AsStringAndSize;
    auto py_as_bytes = PyString_FromStringAndSize;
#endif
}

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

    /* return a list of all the plugins */
    void register_all(object& func) {
        MOLFILE_REGISTER_ALL(&func, register_cb); 
    }

    void destructor(PyObject* obj) { Py_XDECREF(obj); }
    typedef std::shared_ptr<PyObject> objptr;

    DtrWriter* dtr_init(std::string const& path,
                        uint32_t natoms,
                        int mode,
                        uint32_t fpf) {
        return new DtrWriter(path, DtrWriter::Type::DTR, natoms, DtrWriter::Mode(mode), fpf);
    }

    void convert_keyvals_to_keymap(dict keyvals, dtr::KeyMap& map) {
        list items = keyvals.items();
        for (unsigned i=0, n=len(items); i<n; i++) {
            object item = items[i];
            std::string key = extract<std::string>(item[0]);
            object val = item[1];
            PyObject* ptr = val.ptr();
            dtr::Key keyval;
            if (py_string_check(ptr)) {
                const char* s = py_as_string(ptr);
                keyval.set(py_as_string(ptr), strlen(s));
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
    }

    void dtr_append(DtrWriter& w, double time, dict keyvals) {
        dtr::KeyMap keymap;
        convert_keyvals_to_keymap(keyvals, keymap);
        w.append(time, keymap);
    }

    PyObject* py_frame_as_bytes(dict keyvals, bool use_padding) {
        dtr::KeyMap keymap;
        convert_keyvals_to_keymap(keyvals, keymap);
        void* buf = nullptr;
        size_t len = dtr::ConstructFrame(keymap, &buf, use_padding);
        char* ptr = static_cast<char *>(buf);
        PyObject* bytes = py_as_bytes(ptr, len);
        free(buf);
        return bytes;
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
                     arg("frames_per_file")=0)))
            .def("append", dtr_append,
                    (arg("time"),
                     arg("keyvals")))
            .def("sync", &DtrWriter::sync)
            ;

    }

}

BOOST_PYTHON_MODULE(_molfile) {

    _import_array();
    if (PyErr_Occurred()) return;

    docstring_options doc_options;
    doc_options.enable_user_defined();
    doc_options.enable_signatures();
    doc_options.disable_cpp_signatures();

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
                path.c_str(), jiffies_to_ps(key.jiffies()), key.offset(), key.size(),
                first, last+1, filesize, dtrpath.c_str(), dtrsize );
    }

    void py_keyvals(dtr::KeyMap const& keymap, PyObject *dict) {
        for (auto it=keymap.begin(), e=keymap.end(); it!=e; ++it) {
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
                    arr = py_string_from_string_size(
                        reinterpret_cast<const char*>(val.data), val.count);
                    break;
                case dtr::Key::TYPE_UCHAR:
                    arr = py_bytes_from_string_size(
                        reinterpret_cast<const char*>(val.data), val.count);
                    break;
                default:;
            }
            if (!arr) continue;
            char* key = const_cast<char*>(it->first.c_str());
            int rc = PyMapping_SetItemString(dict, key, arr);
            Py_DECREF(arr);
            if (rc==-1) throw_error_already_set();
        }
    }

    dict py_frame_from_bytes(PyObject* bufferobj) {
        Py_buffer view[1];
        if (PyObject_GetBuffer(bufferobj, view, PyBUF_ND)) {
            throw error_already_set();
        }
        std::shared_ptr<Py_buffer> ptr(view, PyBuffer_Release);
        bool swap_endian = false;   // FIXME ignored
        auto keymap = dtr::ParseFrame(view->len, view->buf, &swap_endian);
        dict d;
        py_keyvals(keymap, d.ptr());
        return d;
    }

    const char* keyvals_doc =
        "keyvals(index) -> dict()\n"
        "Read raw fields from frame.\n";

    dict wrap_keyvals(FrameSetReader& self, Py_ssize_t index) {
       const DtrReader *comp = self.component(index);
        if (!comp) {
            PyErr_SetString(PyExc_IndexError, "index out of bounds");
            throw error_already_set();
        }
        void* keybuf = NULL;
        dtr::KeyMap keymap = comp->frame(index, NULL, &keybuf);
        std::shared_ptr<void> dtor(keybuf, free);
        dict d;
        py_keyvals(keymap, d.ptr());
        return d;
    }

    const char * frame_doc = 
        "frame(index, bytes=None, keyvals=None) -> Frame\n"
        "Read bytes from disk if bytes are not provided\n"
        "If keyvals is not None, it should be a dict, and raw data from\n"
        "the frame will be provided.\n";

    object frame(FrameSetReader& self, Py_ssize_t index, object& bytes,
            object keyvals) {
        /* The component method modifies index.  What a dumb API. */
        Py_ssize_t global_index = index;
        const DtrReader *comp = self.component(index);
        if (!comp) {
            PyErr_SetString(PyExc_IndexError, "index out of bounds");
            throw error_already_set();
        }
        Frame *frame = new Frame(comp->natoms(), self.has_velocities(), false);
        dtr::KeyMap keymap;
        void* keybuf = NULL;
        void** keybufptr = keyvals.ptr() == Py_None ? NULL : &keybuf;
        std::shared_ptr<void*> keybuf_dtor(keybufptr, [](void **v) { if (v) free(*v); });
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
            if (py_as_string_size(bytes.ptr(), &data, &size)) {
                delete frame;
                throw error_already_set();
            }
            try {
                comp->frame_from_bytes(data, size, *frame);
                frame->setTime(jiffies_to_ps(comp->keys[index].jiffies()));
            }
            catch (std::exception &e) {
                delete frame;
                PyErr_Format(PyExc_RuntimeError, "Failed parsing frame data of size %ld: %s", size, e.what());
                throw_error_already_set();
            }
        }
        py_keyvals(keymap, keyvals.ptr());

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
        double operator[](boost::python::ssize_t i) const {
            double x;
            self.times(i,1,&x);
            return x;
        }
    };

    boost::python::ssize_t index_near( const FrameSetReader& self, double T ) {
        return ff::at_time_near(self.size(), Oracle(self), T);
    }
    boost::python::ssize_t index_le( const FrameSetReader& self, double T ) {
        return ff::at_time_le(self.size(), Oracle(self), T);
    }
    boost::python::ssize_t index_lt( const FrameSetReader& self, double T ) {
        return ff::at_time_lt(self.size(), Oracle(self), T);
    }
    boost::python::ssize_t index_ge( const FrameSetReader& self, double T ) {
        return ff::at_time_ge(self.size(), Oracle(self), T);
    }
    boost::python::ssize_t index_gt( const FrameSetReader& self, double T ) {
        return ff::at_time_gt(self.size(), Oracle(self), T);
    }
    boost::python::ssize_t total_bytes( const FrameSetReader& self ) {
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
        std::shared_ptr<PyThreadState> rst(_save, PyEval_RestoreThread);
        tk.init(path);
    }

    double tk_interval(Timekeys const& tk) {
        return jiffies_to_ps(tk.interval_jiffies());
    }

}

void desres::molfile::export_dtrreader() {

    class_<Timekeys>("Timekeys", init<>())
        .def("init", tk_init)
        .def("size",                    &Timekeys::size)
        .add_property("framesperfile",  &Timekeys::framesperfile)
        .add_property("framesize",      &Timekeys::framesize)
        .add_property("interval",       tk_interval)
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
                ,arg("keyvals")=object()))
        .def("keyvals", wrap_keyvals, keyvals_doc)
        .def("reload", reload, reload_doc)
        .def("times", get_times)
        ;

    def("dtr_frame_from_bytes", py_frame_from_bytes);
    def("dtr_frame_as_bytes", py_frame_as_bytes,
            (arg("keyvals"), arg("use_padding")=false));
                
    scope().attr("dtr_serialized_version")=dtr_serialized_version();
}

