#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include "graph.hxx"

using namespace desres::msys;
using namespace pybind11;

namespace desres { namespace msys { 

    void export_graph(module m) {

        class_<Graph, GraphPtr>(m, "GraphPtr")
            .def("__eq__", [](Graph* self, Graph* other) { return self==other; })
            .def("__ne__", [](Graph* self, Graph* other) { return self!=other; })
            .def("__hash__", [](Graph* g) { return size_t(g); })
            .def_static("create",  [](SystemPtr mol, IdList const& atoms) { return Graph::create(mol, atoms); })
            .def_static("create_with_colors",  [](SystemPtr mol, IdList const& atoms, IdList const& colors) { return Graph::create(mol, atoms, colors); })
            .def("hash", [](Graph const& g) { return Graph::hash(g.system(), g.atoms()); })
            .def_static("hash_atoms", &Graph::hash)
            .def("size", &Graph::size)
            .def("atoms", &Graph::atoms)
            .def("system", &Graph::system)
            .def("match", [](GraphPtr self, GraphPtr other) { Graph::MatchList m; self->match(other, m); return m; })
            .def("matchAll", [](GraphPtr self, GraphPtr other, bool sub) { std::vector<Graph::MatchList> m; self->matchAll(other, m, sub); return m; })
            ;
    }
}}

