#include "molfilemodule.hxx"
#include "dtrplugin.hxx"
#include "findframe.hxx"

#include <boost/python.hpp>
#include <sstream>

using namespace desres::molfile;
namespace bp=boost::python;
namespace ff=findframe;
using bp::object;
using bp::handle;
using bp::extract;
using bp::error_already_set;
using bp::tuple;
using bp::make_tuple;
using bp::class_;
using bp::no_init;
using bp::arg;
using bp::scope;
using bp::default_call_policies;

namespace {

    FrameSetReader * dtr_from_path( object& pathobj, bool sequential ) {
        std::string path = extract<std::string>(pathobj);
        FrameSetReader * reader = NULL;
        if (StkReader::recognizes(path)) {
            reader = new StkReader(path);
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

    FrameSetReader * dtr_from_bytes( object& byteobj ) {
        const char *buffer;
        Py_ssize_t buffer_len;
        if (PyObject_AsCharBuffer(byteobj.ptr(), &buffer, &buffer_len)) {
            throw error_already_set();
        }
        std::string str(buffer, buffer+buffer_len); // null-terminate 
        /* check the path to see if we have an stk or not */
        const char * space = strchr(str.c_str(), ' ');
        if (!space || space-str.c_str() < 4) {
            PyErr_SetString(PyExc_ValueError, "misformatted bytes");
            throw error_already_set();
        }
        std::string extension(space-4, space);
        FrameSetReader * reader = NULL;
        if (extension==".stk") {
            reader = new StkReader("<from bytes>");
        } else {
            reader = new DtrReader("<from bytes>");
        }
        std::istringstream in(str);
        reader->load(in);
        if (!in) {
            delete reader;
            PyErr_SetString(PyExc_ValueError, "reading bytes failed");
            throw error_already_set();
        }
        return reader;
    }

    FrameSetReader * dtr_init( object& path, object& bytes, bool sequential ) {
        if (!path.is_none()) return dtr_from_path(path, sequential);
        else if (!bytes.is_none()) return dtr_from_bytes(bytes);
        else {
            PyErr_Format(PyExc_ValueError, "Must supply either path or bytes");
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
        "frame(index, bytes=None, with_gids=False) -> Frame object\n"
        "Read bytes from disk if bytes are not provided\n";

    object frame(const FrameSetReader& self, Py_ssize_t index, object& bytes,
            bool with_gids) {
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
        if (bytes.is_none()) {
            PyThreadState *_save;
            _save = PyEval_SaveThread();
            try {
                comp->frame(index, *frame);
            }
            catch (std::exception &e) {
                PyEval_RestoreThread(_save);
                delete frame;
                PyErr_Format(PyExc_IOError, "Error reading frame: global index %ld dtr path %s local index %ld frame file %s\n%s",
                        global_index, comp->path().c_str(),
                        index, comp->framefile(index).c_str(),
                        e.what());
                bp::throw_error_already_set();
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
                bp::throw_error_already_set();
            }
        }
        return object(frame);
    }

    const char * dump_doc = 
        "dump() -> serialized version of self\n";

    object dump(const FrameSetReader& self) {
        std::ostringstream out;
        self.dump(out);
        const std::string &str = out.str();
        return object(handle<>(
                    PyString_FromStringAndSize(str.c_str(), str.size())));
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

    std::string my_path(FrameSetReader const& self) {
        return self.path();
    }

    FrameSetReader* from_timekeys(std::string const& path,
                                  bp::list lfnames, bp::list ltimekeys) {
        if (bp::len(lfnames) != bp::len(ltimekeys)) {
            PyErr_Format(PyExc_ValueError, "fnames and timekeys must be the same length");
            bp::throw_error_already_set();
        }
        unsigned i,n = bp::len(lfnames);
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

    class_<Timekeys>("Timekeys", bp::init<>())
        .def("init", tk_init)
        ;

    class_<FrameSetReader, boost::noncopyable>("DtrReader", no_init)
        .def("__init__", 
                make_constructor( 
                    dtr_init,
                    default_call_policies(),
                    (arg("path")=object(), 
                     arg("bytes")=object(),
                     /* WARNING: use sequential=True only from one thread! */
                     arg("sequential")=false )))
        .add_property("path", my_path)
        .add_property("natoms", &FrameSetReader::natoms)
        .add_property("nframes", &FrameSetReader::size)
        .add_property("nframesets", &FrameSetReader::nframesets)

        .def("fromTimekeys", from_timekeys, 
                bp::return_value_policy<bp::manage_new_object>())
        .staticmethod("fromTimekeys")

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
                ,arg("with_gids")=false))
        .def("dump", dump, dump_doc)
        .def("reload", reload, reload_doc)
        .def("times", get_times)
        ;

                
    scope().attr("dtr_serialized_version")=dtr_serialized_version();
}

