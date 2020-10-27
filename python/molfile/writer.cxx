#include "molfilemodule.hxx"
#include <numpy/arrayobject.h>

namespace {
    using namespace desres::molfile;
    void write_grid(Writer& w, dict d, object data) {
        grid_t g;
        memset(&g,0,sizeof(g));

        /* dataname */
        strncpy(g.dataname, d["name"].cast<std::string>().data(), sizeof(g.dataname));
        g.dataname[sizeof(g.dataname)-1]='\0';

        /* origin */
        list origin = d["origin"];
        for (int i=0; i<3; i++) g.origin[i] = origin[i].cast<float>();

        /* axis */
        list xaxis = d["xaxis"];
        list yaxis = d["yaxis"];
        list zaxis = d["zaxis"];
        for (int i=0; i<3; i++) {
            g.xaxis[i] = xaxis[i].cast<float>();
            g.yaxis[i] = yaxis[i].cast<float>();
            g.zaxis[i] = zaxis[i].cast<float>();
        }

        /* dims -> size */
        list size = d["dims"];
        g.xsize = size[0].cast<int>();
        g.ysize = size[1].cast<int>();
        g.zsize = size[2].cast<int>();

        w.write_grid(g, (const float *)PyArray_DATA(data.ptr()));
    }
}

void desres::molfile::export_writer(module m) {

    class_<Writer>(m, "Writer", "Structure or trajectory open for writing")
        .def_property_readonly("natoms", &Writer::natoms)
        .def_property_readonly("path", &Writer::path)
        .def_property_readonly("can_sync", &Writer::can_sync)
        .def("sync", &Writer::sync)
        .def("truncate", &Writer::truncate)
        .def("close", &Writer::close)
        .def("frame", &Writer::write_frame, return_value_policy::reference)
        .def("_grid", write_grid)
        ;
}

