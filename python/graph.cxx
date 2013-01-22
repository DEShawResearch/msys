#include "graph.hxx"
#include "wrap_obj.hxx"

using namespace desres::msys;

namespace {
    object match(Graph const& self, GraphPtr other) {
        std::vector<IdPair> matches;
        if (self.match(other, matches)) {
            list L;
            BOOST_FOREACH(IdPair const& p, matches) {
                L.append(make_tuple(p.first, p.second));
            }
            return L;
        }
        return object();
    }

    object matchAll(Graph const& self, GraphPtr other) {
        std::vector<std::vector<IdPair> > matches;
        self.matchAll(other, matches);
        list outer_L;
        BOOST_FOREACH(std::vector<IdPair> const& v, matches) {
            list L;
            BOOST_FOREACH(IdPair const& p, v)
                L.append(make_tuple(p.first, p.second));
            outer_L.append(L);
        }
        return outer_L;
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
            .def("atoms", &Graph::atoms, return_const())
            .def("match", match)
            .def("matchAll", matchAll)
            ;
    }
}}
