#include "molfilemodule.hxx"
#include <pybind11/stl.h>

using namespace desres::molfile;

namespace {

    handle frame_pos(object& obj) {
        auto& self = obj.cast<Frame&>();
        if (!self.pos()) return none();
        Py_ssize_t dims[2] = { (Py_ssize_t)self.natoms(), 3 };
        return ((
                    backed_vector(2, dims, desres::molfile::FLOAT, self.pos(), obj.ptr())));
    }

    handle frame_vel(object& obj) {
        auto& self = obj.cast<Frame&>();
        if (!self.vel()) return none();
        Py_ssize_t dims[2] = { (Py_ssize_t)self.natoms(), 3 };
        return (( 
                    backed_vector(2, dims, desres::molfile::FLOAT, self.vel(), obj.ptr())));
    }

    handle frame_dpos(object& obj) {
        auto& self = obj.cast<Frame&>();
        if (!self.dpos()) return none();
        Py_ssize_t dims[2] = { (Py_ssize_t)self.natoms(), 3 };
        return ((
                    backed_vector(2, dims, desres::molfile::DOUBLE, self.dpos(), obj.ptr())));
    }

    handle frame_dvel(object& obj) {
        auto& self = obj.cast<Frame&>();
        if (!self.dvel()) return none();
        Py_ssize_t dims[2] = { (Py_ssize_t)self.natoms(), 3 };
        return (( 
                    backed_vector(2, dims, desres::molfile::DOUBLE, self.dvel(), obj.ptr())));
    }

    handle frame_box(object& obj) {
        auto& self = obj.cast<Frame&>();
        Py_ssize_t dims[2] = {3,3};
        return (( 
                    backed_vector(2, dims, desres::molfile::DOUBLE, self.box(), obj.ptr())));
    }

    handle frame_ptensor(object& obj) {
        auto& self = obj.cast<Frame&>();
        Py_ssize_t dims[2] = {3,3};
        return (( 
                    backed_vector(2, dims, desres::molfile::DOUBLE, self.pressure_tensor(), 
                      obj.ptr())));
    }

    handle frame_virial(object& obj) {
        auto& self = obj.cast<Frame&>();
        Py_ssize_t dims[2] = {3,3};
        return (( 
                    backed_vector(2, dims, desres::molfile::DOUBLE, self.virial_tensor(), 
                      obj.ptr())));
    }

    Frame * frame_select( const Frame& self, std::vector<unsigned> const& indices) {
        Frame * frame = new Frame(indices.size(), self.vel()!=NULL);
        float * pos = frame->pos();
        float * vel = frame->vel();
        try {
            for (auto ind : indices) {
                if (ind<0 || ind>=(Py_ssize_t)self.natoms()) {
                    PyErr_Format(PyExc_ValueError, 
                            "Index %ld is out of range", ind);
                    throw error_already_set();
                }
                for (int j=0; j<3; j++)          *pos++=self.pos()[3*ind+j];
                if (vel) for (int j=0; j<3; j++) *vel++=self.vel()[3*ind+j];
            }
        } catch (std::exception& e) {
            delete frame;
            throw;
        }
        return frame;
    }

    void frame_moveby( Frame& self, float x, float y, float z ) {
        float *p=self.pos();
        const float * const end = p+3*self.natoms();
        while (p!=end) {
            p[0] += x;
            p[1] += y;
            p[2] += z;
            p += 3;
        }
    }

#define WRAP_SCALAR_GET( x ) \
    object get_##x ( const Frame& f ) { \
        if (f.has_##x ()) return cast(f. x ()) ; \
        return none(); \
    }

    WRAP_SCALAR_GET(total_energy)
    WRAP_SCALAR_GET(potential_energy)
    WRAP_SCALAR_GET(kinetic_energy)
    WRAP_SCALAR_GET(extended_energy)
    WRAP_SCALAR_GET(temperature)
    WRAP_SCALAR_GET(pressure)
}

void desres::molfile::export_frame(module m) {

#define SCALARPROP(x) .def_property(#x, get_##x, &Frame:: set_##x)

    class_<Frame>(m, "Frame")
        .def(init<size_t, bool, bool>(), 
                     arg("natoms")
                    ,arg("with_velocities")=false
                    ,arg("double_precision")=false)
        .def_property("time", &Frame::time, &Frame::setTime)
        .def_property_readonly("position",   frame_pos)
        .def_property_readonly("pos",        frame_pos)
        .def_property_readonly("fpos",       frame_pos)
        .def_property_readonly("dpos",       frame_dpos)
        .def_property_readonly("velocity",   frame_vel)
        .def_property_readonly("vel",        frame_vel)
        .def_property_readonly("fvel",       frame_vel)
        .def_property_readonly("dvel",       frame_dvel)
        .def_property_readonly("box", frame_box)
        .def_property_readonly("pressure_tensor", frame_ptensor)
        .def_property_readonly("virial_tensor", frame_virial)
        SCALARPROP(total_energy)
        SCALARPROP(potential_energy)
        SCALARPROP(kinetic_energy)
        SCALARPROP(extended_energy)
        SCALARPROP(temperature)
        SCALARPROP(pressure)
        .def("select", frame_select, arg("indices"))
        .def("moveby", frame_moveby, arg("x"), arg("y"), arg("z"))
        ;
}

