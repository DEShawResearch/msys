#include "wrap_obj.hxx"

namespace desres { namespace msys { 

    void export_vector() {

        class_<Vec3>("Vec3", init<Float,Float,Float>())
            .add_property("x", &Vec3::x, &Vec3::x)
            .add_property("y", &Vec3::y, &Vec3::y)
            .add_property("z", &Vec3::z, &Vec3::z)
            .def("__eq__", list_eq<Vec3>)
            .def("__ne__", list_ne<Vec3>)
            ;

        declare_list<IdList>("IdList");
        declare_list<std::vector<Int> >("IntVec");
        declare_list<std::vector<Float> >("FloatVec");
        declare_list<std::vector<String> >("StringVec");

        def("bad", bad);
        scope().attr("BadId") = (Id)BadId;
    }

}}
