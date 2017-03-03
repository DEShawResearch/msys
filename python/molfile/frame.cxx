#include "molfilemodule.hxx"
#include <boost/python.hpp>

using namespace desres::molfile;
using namespace boost::python;

namespace {

    object frame_pos(object& obj) {
        Frame& self = extract<Frame&>(obj);
        if (!self.pos()) return object();
        Py_ssize_t dims[2] = { (Py_ssize_t)self.natoms(), 3 };
        return object(handle<>(
                    backed_vector(2, dims, desres::molfile::FLOAT, self.pos(), obj.ptr())));
    }

    object frame_vel(object& obj) {
        Frame& self = extract<Frame&>(obj);
        if (!self.vel()) return object();
        Py_ssize_t dims[2] = { (Py_ssize_t)self.natoms(), 3 };
        return object(handle<>( 
                    backed_vector(2, dims, desres::molfile::FLOAT, self.vel(), obj.ptr())));
    }

    object frame_dpos(object& obj) {
        Frame& self = extract<Frame&>(obj);
        if (!self.dpos()) return object();
        Py_ssize_t dims[2] = { (Py_ssize_t)self.natoms(), 3 };
        return object(handle<>(
                    backed_vector(2, dims, desres::molfile::DOUBLE, self.dpos(), obj.ptr())));
    }

    object frame_dvel(object& obj) {
        Frame& self = extract<Frame&>(obj);
        if (!self.dvel()) return object();
        Py_ssize_t dims[2] = { (Py_ssize_t)self.natoms(), 3 };
        return object(handle<>( 
                    backed_vector(2, dims, desres::molfile::DOUBLE, self.dvel(), obj.ptr())));
    }

    object frame_box(object& obj) {
        Frame& self = extract<Frame&>(obj);
        Py_ssize_t dims[2] = {3,3};
        return object(handle<>( 
                    backed_vector(2, dims, desres::molfile::DOUBLE, self.box(), obj.ptr())));
    }

    object frame_gid(object& obj) {
        Frame& self = extract<Frame&>(obj);
        if (!self.gid()) return object();
        Py_ssize_t dims[1] = {(Py_ssize_t)self.natoms() };
        return object(handle<>( 
                    backed_vector(1, dims, desres::molfile::INT, self.gid(), obj.ptr())));
    }

    object frame_ptensor(object& obj) {
        Frame& self = extract<Frame&>(obj);
        Py_ssize_t dims[2] = {3,3};
        return object(handle<>( 
                    backed_vector(2, dims, desres::molfile::DOUBLE, self.pressure_tensor(), 
                      obj.ptr())));
    }

    object frame_virial(object& obj) {
        Frame& self = extract<Frame&>(obj);
        Py_ssize_t dims[2] = {3,3};
        return object(handle<>( 
                    backed_vector(2, dims, desres::molfile::DOUBLE, self.virial_tensor(), 
                      obj.ptr())));
    }

    Frame * frame_select( const Frame& self, object& indices ) {
        long i,n = len(indices);
        Frame * frame = new Frame(n, self.vel()!=NULL);
        float * pos = frame->pos();
        float * vel = frame->vel();
        try {
            for (i=0; i<n; i++) {
                Py_ssize_t ind = extract<Py_ssize_t>(indices[i]);
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
        if (f.has_##x ()) return object(f. x ()) ; \
        return object(); \
    }

    WRAP_SCALAR_GET(total_energy)
    WRAP_SCALAR_GET(potential_energy)
    WRAP_SCALAR_GET(kinetic_energy)
    WRAP_SCALAR_GET(extended_energy)
    WRAP_SCALAR_GET(temperature)
    WRAP_SCALAR_GET(pressure)
}

static Frame* new_frame(size_t n, bool vel, bool dbl, bool with_gid) {
    return new Frame(n,vel,dbl,with_gid);
}

void desres::molfile::export_frame() {

#define SCALARPROP(x) .add_property(#x, get_##x, &Frame:: set_##x)

    class_<Frame>("Frame", no_init)
        .def("__init__", make_constructor(
                    new_frame,
                    default_call_policies(),
                    (arg("natoms")
                    ,arg("with_velocities")=false
                    ,arg("double_precision")=false
                    ,arg("with_gids")=false
                    )))
        .add_property("time", &Frame::time, &Frame::setTime)
        .add_property("position",   frame_pos)
        .add_property("pos",        frame_pos)
        .add_property("fpos",       frame_pos)
        .add_property("dpos",       frame_dpos)
        .add_property("velocity",   frame_vel)
        .add_property("vel",        frame_vel)
        .add_property("fvel",       frame_vel)
        .add_property("dvel",       frame_dvel)
        .add_property("box", frame_box)
        .add_property("gid", frame_gid)
        .add_property("pressure_tensor", frame_ptensor)
        .add_property("virial_tensor", frame_virial)
        SCALARPROP(total_energy)
        SCALARPROP(potential_energy)
        SCALARPROP(kinetic_energy)
        SCALARPROP(extended_energy)
        SCALARPROP(temperature)
        SCALARPROP(pressure)
        .def("select", frame_select, 
                (arg("indices")),
                return_value_policy<manage_new_object>())
        .def("moveby", frame_moveby,
                (arg("x"), arg("y"), arg("z")))
        ;
}

