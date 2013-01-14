#include "graph.hxx"
#include "wrap_obj.hxx"

using namespace desres::msys;

namespace {
    inline bool sort_by_second(IdPair const& a, IdPair const& b) { 
        return a.second<b.second; 
    }
    object match(Graph const& self, GraphPtr other) {
        std::vector<IdPair> matches;
        if (self.match(other, matches)) {
            std::sort(matches.begin(), matches.end(), sort_by_second);
            list L;
            BOOST_FOREACH(IdPair const& p, matches) {
                L.append(p.first);
            }
            return L;
        }
        return object();
    }
}

namespace desres { namespace msys { 

    void export_graph() {

        class_<Graph, GraphPtr, boost::noncopyable>("GraphPtr", no_init)
            .def("__eq__",      list_eq<GraphPtr>)
            .def("__ne__",      list_ne<GraphPtr>)
            .def("__hash__",    obj_hash<GraphPtr>)
            .def("create",  &Graph::create).staticmethod("create")
            .def("hash", &Graph::hash)
            .def("size", &Graph::size)
            .def("match", match)
            ;
    }
}}

