#include "wrap_obj.hxx"

namespace desres { namespace msys { 

    void export_vector() {

        class_<IdList>("IdList")
            .def(vector_indexing_suite<IdList>())
            /* for some reason, operator== returns false even when all
             * the component elements compare equal.  What's that about? */
            .def("__eq__", list_eq<IdList>)
            .def("__ne__", list_ne<IdList>)
            ;

        class_<MultiIdList>("MultiIdList")
            .def(vector_indexing_suite<MultiIdList>())
            /* for some reason, operator== returns false even when all
             * the component elements compare equal.  What's that about? */
            .def("__eq__", list_eq<MultiIdList>)
            .def("__ne__", list_ne<MultiIdList>)
            ;
    
    
        def("bad", bad);
        scope().attr("BadId") = (Id)BadId;
    }

}}
