#include <pybind11/pybind11.h>
#include "molfilemodule.hxx"
#include "numpy/arrayobject.h"
#include <map>
#include <string>
#include <sstream>

#include <msys/molfile/findframe.hxx>
#include <msys/molfile/dtrplugin.hxx>
#include <msys/molfile/dtrframe.hxx>

using namespace desres::molfile;
using namespace pybind11;
namespace ff=findframe;

namespace {
    auto py_string_check = [](PyObject* o) { return PyUnicode_Check(o); };
    auto py_string_from_string_size = PyUnicode_FromStringAndSize;
    auto py_as_string = PyUnicode_AsUTF8;
    auto py_as_string_size = PyBytes_AsStringAndSize;
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

    void convert_keyvals_to_keymap(dict keyvals, dtr::KeyMap& map) {
        for (auto item : keyvals) {
            auto key = item.first.cast<std::string>();
            PyObject* ptr = item.second.ptr();
            dtr::Key keyval;
            if (py_string_check(ptr)) {
                const char* s = py_as_string(ptr);
                keyval.set(py_as_string(ptr), strlen(s));
            } else if (PyByteArray_Check(ptr)) {
                keyval.set((unsigned char *)PyByteArray_AsString(ptr), PyByteArray_Size(ptr));
            } else if (PyArray_Check(ptr)) {
                if (PyArray_NDIM(ptr)!=1) {
                    PyErr_Format(PyExc_ValueError, "array for key %s must be 1d", key.c_str());
                    throw error_already_set();
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
                        throw error_already_set();
                }
            } else {
                PyErr_Format(PyExc_ValueError, "values must be string, bytearray, or numpy array");
                throw error_already_set();
            }
            map[key] = keyval;
        }
    }

    bytes py_frame_as_bytes(dict keyvals, bool use_padding, double precision) {
        dtr::KeyMap keymap;
        convert_keyvals_to_keymap(keyvals, keymap);
        void* buf = nullptr;
        size_t len = dtr::ConstructFrame(keymap, &buf, use_padding, precision);
        char* ptr = static_cast<char *>(buf);
        auto obj = bytes(ptr, len);
        free(buf);
        return obj;
    }

}

namespace {

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
                case dtr::Key::TYPE_UCHAR:
                    arr = PyByteArray_FromStringAndSize((const char*)val.data, val.count);
                    break;
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
                    if (!arr) {
                        PyErr_Clear();
                        arr = PyByteArray_FromStringAndSize((const char*)val.data, val.count);
                    }
                    break;
                default:;
            }
            if (!arr) continue;
            char* key = const_cast<char*>(it->first.c_str());
            int rc = PyMapping_SetItemString(dict, key, arr);
            Py_DECREF(arr);
            if (rc==-1) throw error_already_set();
        }
    }

    dict py_frame_from_bytes(buffer buf) {
        auto req = buf.request();
        bool swap_endian = false;   // FIXME ignored
        void* allocated = nullptr;  // allocated when buffer contains compressed positions
        auto keymap = dtr::ParseFrame(req.size, req.ptr, &swap_endian, &allocated);
        dict d;
        py_keyvals(keymap, d.ptr());
        free(allocated);    // py_keyvals makes a copy
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

    std::unique_ptr<Frame> frame(FrameSetReader& self, Py_ssize_t index, object& bytes,
            object keyvals) {
        /* The component method modifies index.  What a dumb API. */
        Py_ssize_t global_index = index;
        const DtrReader *comp = self.component(index);
        if (!comp) {
            PyErr_SetString(PyExc_IndexError, "index out of bounds");
            throw error_already_set();
        }
        auto frame = std::unique_ptr<Frame>(new Frame(comp->natoms(), self.has_velocities(), false));
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
                PyErr_Format(PyExc_IOError, "Error reading frame: global index %ld dtr path %s local index %ld frame file %s\n%s",
                        global_index, comp->path().c_str(),
                        index, comp->framefile(index).c_str(),
                        e.what());
                throw error_already_set();
            }
            PyEval_RestoreThread(_save);

        } else {
            Py_ssize_t size;
            char * data;
            if (py_as_string_size(bytes.ptr(), &data, &size)) {
                throw error_already_set();
            }
            try {
                comp->frame_from_bytes(data, size, *frame);
                frame->setTime(jiffies_to_ps(comp->keys[index].jiffies()));
            }
            catch (std::exception &e) {
                PyErr_Format(PyExc_RuntimeError, "Failed parsing frame data of size %ld: %s", size, e.what());
                throw error_already_set();
            }
        }
        py_keyvals(keymap, keyvals.ptr());

        return frame;
    }

    const char reload_doc[] =
        "reload() -> number of timekeys reloaded -- reload frames in the dtr/stk";
    int reload(FrameSetReader& self) {
        int changed=0;
        self.init(&changed);
        return changed;
    }

    handle get_times(const FrameSetReader& self) {
        Py_ssize_t n = self.size();
        Py_ssize_t dims[1] = { n };
        PyObject * arr = backed_vector( 1, dims, DOUBLE, NULL, NULL );
        if (self.times(0,n,(double *)array_data(arr))!=n) {
            Py_DECREF(arr);
            arr=NULL;
            PyErr_Format(PyExc_RuntimeError, "Error reading times");
            throw error_already_set();
        }
        return handle(arr);
    }
    ssize_t frameset_size(const FrameSetReader& self, Py_ssize_t n) {
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

    std::string my_path(FrameSetReader const& self) {
        return self.path();
    }

    std::unique_ptr<FrameSetReader> from_timekeys(std::string const& path,
            std::vector<std::string> fnames, std::vector<Timekeys> timekeys) {
        if (fnames.size() != timekeys.size()) {
            PyErr_Format(PyExc_ValueError, "fnames and timekeys must be the same length");
            throw error_already_set();
        }
        std::unique_ptr<StkReader> stk(new StkReader(path));
        stk->append(fnames, timekeys);
        return std::move(stk);
    }

    double tk_interval(Timekeys const& tk) {
        return jiffies_to_ps(tk.interval_jiffies());
    }

    dict frameset_metadata(FrameSetReader& r) {
        dict d;
        if (r.nframesets()>0) {
            auto meta = r.frameset(0)->get_meta();
            py_keyvals(*meta->get_frame_map(), d.ptr());
        }
        return d;
    }

}

PYBIND11_MODULE(_molfile, m) {
    _import_array();
    if (PyErr_Occurred()) throw error_already_set();

    export_plugin(m);
    export_reader(m);
    export_frame(m);
    export_writer(m);

    // I could use a static to initialize these types on first access,
    // but then the type objects wouldn't be present in the module
    // until first use.  This would cause help(molfile) to be incomplete.
    if (initialize_atom()) return;

    static PyObject *pyAtom = reinterpret_cast<PyObject *>(&AtomType);
    Py_INCREF(pyAtom);
    m.attr("Atom")=handle(pyAtom);

    MOLFILE_INIT_ALL;
    m.def("register_all", register_all);

    auto writer = class_<DtrWriter>(m, "DtrWriter");
    enum_<DtrWriter::Type>(writer, "Type")
        .value("DTR", DtrWriter::Type::DTR)
        .value("ETR", DtrWriter::Type::ETR)
        .export_values()
        ;

    writer
        .def(init([](std::string const& path, uint32_t natoms, int mode, uint32_t fpf,
                     DtrWriter::Type type, double precision, object metadata) {
                std::unique_ptr<dtr::KeyMap> keymap;
                if (!metadata.is_none()) {
                    keymap.reset(new dtr::KeyMap);
                    convert_keyvals_to_keymap(dict(metadata), *keymap);
                }
                return new DtrWriter(path, type, natoms, DtrWriter::Mode(mode), fpf, keymap.get(), precision);
            }), arg("path"), arg("natoms"), arg("mode")=0, arg("frames_per_file")=0,
                arg("format")=DtrWriter::Type::DTR, arg("precision")=0.0, arg("metadata")=none())
        .def("append", [](DtrWriter& w, double time, dict keyvals) {
            dtr::KeyMap keymap;
            convert_keyvals_to_keymap(keyvals, keymap);
            w.append(time, keymap);
            }, arg("time"), arg("keyvals"))
        .def("sync", &DtrWriter::sync)
        .def("close", &DtrWriter::close)
        ;


    class_<Timekeys>(m, "Timekeys")
        .def(init<>())
        .def("init", &Timekeys::init)
        .def("size", &Timekeys::size)
        .def_property_readonly("framesperfile",  &Timekeys::framesperfile)
        .def_property_readonly("framesize",      &Timekeys::framesize)
        .def_property_readonly("interval",       tk_interval)
        ;

    class_<FrameSetReader>(m, "DtrReader")
        .def(init([](std::string const& path, bool sequential) {
            std::unique_ptr<FrameSetReader> r;
            if (StkReader::recognizes(path)) {
                r.reset(new StkReader(path, sequential ? DtrReader::SequentialAccess : 0));
            } else {
                r.reset(new DtrReader(path, sequential ? DtrReader::SequentialAccess : 0));
            }
            r->init();
            return r;
        }), arg("path"), arg("sequential")=false)
        .def_property_readonly("path", &FrameSetReader::path)
        .def_property_readonly("natoms", &FrameSetReader::natoms)
        .def_property_readonly("nframes", &FrameSetReader::size)
        .def_property_readonly("nframesets", &FrameSetReader::nframesets)
        .def_property_readonly("metadata", frameset_metadata)

        .def_static("fromTimekeys", from_timekeys)
        .def("total_bytes", &FrameSetReader::total_bytes)
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
                arg("index") 
               ,arg("bytes")=none()
               ,arg("keyvals")=none())
        .def("keyvals", wrap_keyvals, keyvals_doc)
        .def("reload", reload, reload_doc)
        .def("times", get_times)
        ;

    m.def("dtr_frame_from_bytes", py_frame_from_bytes);
    m.def("dtr_frame_as_bytes", py_frame_as_bytes,
            arg("keyvals"), arg("use_padding")=false, arg("precision")=0.0);
                
    m.attr("dtr_serialized_version")=dtr_serialized_version();

}

