#include "wrap_obj.hxx"
#include "graph.hxx"

using namespace desres::msys;

namespace {

    GraphPtr create(SystemPtr mol, list atoms) {
        return Graph::create(mol, ids_from_python(atoms));
    }

    GraphPtr create_with_colors(SystemPtr mol, list atoms, list colors) {
        return Graph::create(mol, ids_from_python(atoms), ids_from_python(colors));
    }

    std::string hash(Graph const& self) {
        return Graph::hash(self.system(), self.atoms());
    }

    object match(Graph const& self, GraphPtr other) {
        std::vector<IdPair> matches;
        if (self.match(other, matches)) {
            list L;
            for (IdPair const& p : matches) {
                L.append(make_tuple(p.first, p.second));
            }
            return L;
        }
        return object();
    }

    object matchAll(Graph const& self, GraphPtr other, bool substructure) {
        std::vector<std::vector<IdPair> > matches;
        self.matchAll(other, matches, substructure);
        list outer_L;
        for (std::vector<IdPair> const& v : matches) {
            list L;
            for (IdPair const& p : v)
                L.append(make_tuple(p.first, p.second));
            outer_L.append(L);
        }
        return outer_L;
    }
    list graph_atoms(Graph const& g) {
        return to_python(g.atoms());
    }
}

namespace desres { namespace msys { 

    void export_graph() {

        class_<Graph, GraphPtr, boost::noncopyable>("GraphPtr", no_init)
            .def("__eq__",      list_eq<GraphPtr>)
            .def("__ne__",      list_ne<GraphPtr>)
            .def("__hash__",    obj_hash<GraphPtr>)
            .def("create",  create).staticmethod("create")
            .def("create_with_colors",  create_with_colors).staticmethod("create_with_colors")
            .def("hash", hash)
            .def("size", &Graph::size)
            .def("atoms", graph_atoms)
            .def("system", &Graph::system)
            .def("match", match)
            .def("matchAll", matchAll)
            ;
    }
}}

