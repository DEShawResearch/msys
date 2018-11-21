#include "molfilemodule.hxx"
#include <boost/python.hpp>
#include <numpy/arrayobject.h>

using namespace boost::python;

namespace {
    using namespace desres::molfile;
    void write_grid(Writer& w, dict d, PyObject* data) {
        grid_t g;
        memset(&g,0,sizeof(g));

        /* dataname */
        strncpy(g.dataname, extract<const char *>(d["name"]), sizeof(g.dataname));
        g.dataname[sizeof(g.dataname)-1]='\0';

        /* origin */
        object origin = d["origin"];
        for (int i=0; i<3; i++) g.origin[i] = extract<float>(origin[i]);

        /* axis */
        object xaxis = d["xaxis"];
        object yaxis = d["yaxis"];
        object zaxis = d["zaxis"];
        for (int i=0; i<3; i++) {
            g.xaxis[i] = extract<float>(xaxis[i]);
            g.yaxis[i] = extract<float>(yaxis[i]);
            g.zaxis[i] = extract<float>(zaxis[i]);
        }

        /* dims -> size */
        object size = d["dims"];
        g.xsize = extract<int>(size[0]);
        g.ysize = extract<int>(size[1]);
        g.zsize = extract<int>(size[2]);

        w.write_grid(g, (const float *)PyArray_DATA(data));
    }
}

void desres::molfile::export_writer() {

    class_<Writer>("Writer", "Structure or trajectory open for writing",no_init)
        .add_property("natoms", &Writer::natoms)
        .add_property("path", &Writer::path)
        .add_property("can_sync", &Writer::can_sync)
        .def("sync", &Writer::sync)
        .def("truncate", &Writer::truncate)
        .def("close", &Writer::close)
        .def("frame", &Writer::write_frame,
                return_internal_reference<1>())
        .def("_grid", write_grid)
        ;
}

