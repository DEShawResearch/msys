#include "sssr.hxx"
#include "dms.hxx"

#include <iostream>

using namespace desres::msys;
using namespace desres::msys::SSSR;

static void test_biconnected_components() {
    GraphRepr graph;
    graph.edges.push_back(Edge(0,1));
    graph.edges.push_back(Edge(0,7));
    graph.edges.push_back(Edge(0,18));
    graph.edges.push_back(Edge(2,3));
    graph.edges.push_back(Edge(2,18));
    graph.edges.push_back(Edge(3,4));
    graph.edges.push_back(Edge(4,5));
    graph.edges.push_back(Edge(4,7));
    graph.edges.push_back(Edge(5,6));
    graph.edges.push_back(Edge(6,17));
    graph.edges.push_back(Edge(7,16));
    graph.edges.push_back(Edge(8,9));
    graph.edges.push_back(Edge(8,17));
    graph.edges.push_back(Edge(9,10));
    graph.edges.push_back(Edge(9,13));
    graph.edges.push_back(Edge(10,11));
    graph.edges.push_back(Edge(10,12));
    graph.edges.push_back(Edge(11,12));
    graph.edges.push_back(Edge(13,14));
    graph.edges.push_back(Edge(14,15));
    graph.edges.push_back(Edge(15,17));
    graph.edges.push_back(Edge(16,17));
    graph.edges.push_back(Edge(19,20));
    graph.v_to_e.resize(21, std::vector<int>());
    graph.v_to_e[0].push_back(0);
    graph.v_to_e[0].push_back(1);
    graph.v_to_e[0].push_back(2);
    graph.v_to_e[1].push_back(0);
    graph.v_to_e[2].push_back(3);
    graph.v_to_e[2].push_back(4);
    graph.v_to_e[3].push_back(3);
    graph.v_to_e[3].push_back(5);
    graph.v_to_e[4].push_back(5);
    graph.v_to_e[4].push_back(6);
    graph.v_to_e[4].push_back(7);
    graph.v_to_e[5].push_back(6);
    graph.v_to_e[5].push_back(8);
    graph.v_to_e[6].push_back(8);
    graph.v_to_e[6].push_back(9);
    graph.v_to_e[7].push_back(1);
    graph.v_to_e[7].push_back(7);
    graph.v_to_e[7].push_back(10);
    graph.v_to_e[8].push_back(11);
    graph.v_to_e[8].push_back(12);
    graph.v_to_e[9].push_back(11);
    graph.v_to_e[9].push_back(13);
    graph.v_to_e[9].push_back(14);
    graph.v_to_e[10].push_back(13);
    graph.v_to_e[10].push_back(15);
    graph.v_to_e[10].push_back(16);
    graph.v_to_e[11].push_back(15);
    graph.v_to_e[11].push_back(17);
    graph.v_to_e[12].push_back(16);
    graph.v_to_e[12].push_back(17);
    graph.v_to_e[13].push_back(14);
    graph.v_to_e[13].push_back(18);
    graph.v_to_e[14].push_back(18);
    graph.v_to_e[14].push_back(19);
    graph.v_to_e[15].push_back(19);
    graph.v_to_e[15].push_back(20);
    graph.v_to_e[16].push_back(10);
    graph.v_to_e[16].push_back(21);
    graph.v_to_e[17].push_back(9);
    graph.v_to_e[17].push_back(12);
    graph.v_to_e[17].push_back(20);
    graph.v_to_e[17].push_back(21);
    graph.v_to_e[18].push_back(2);
    graph.v_to_e[18].push_back(4);
    graph.v_to_e[19].push_back(22);
    graph.v_to_e[20].push_back(22);

    std::cout << "Testing biconnected components" << std::endl;
    std::vector<GraphRepr> components;
    std::vector<std::vector<int> > components_idx;
    get_biconnected_components(graph, components, components_idx);
    for (unsigned i = 0; i < components.size(); ++i) {
        std::cout << "Component " << i << ":" << std::endl;
        for (unsigned j = 0; j < components[i].edges.size(); ++j)
            std::cout << "(" << components_idx[i][components[i].edges[j].first]
                << "," << components_idx[i][components[i].edges[j].second]
                << ") ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

static void test_min_cycle_basis() {
    GraphRepr graph;
    graph.edges.push_back(Edge(0,1));
    graph.edges.push_back(Edge(1,2));
    graph.edges.push_back(Edge(2,3));
    graph.edges.push_back(Edge(3,4));
    graph.edges.push_back(Edge(4,5));
    graph.edges.push_back(Edge(5,0));
    graph.edges.push_back(Edge(0,3));
    graph.edges.push_back(Edge(0,6));
    graph.edges.push_back(Edge(6,1));
    graph.edges.push_back(Edge(1,7));
    graph.edges.push_back(Edge(7,2));
    graph.edges.push_back(Edge(2,8));
    graph.edges.push_back(Edge(8,3));
    graph.edges.push_back(Edge(3,9));
    graph.edges.push_back(Edge(9,4));
    graph.edges.push_back(Edge(4,10));
    graph.edges.push_back(Edge(10,5));
    graph.edges.push_back(Edge(5,11));
    graph.edges.push_back(Edge(11,0));
    graph.edges.push_back(Edge(0,12));
    graph.edges.push_back(Edge(12,3));
    graph.v_to_e.resize(13, std::vector<int>());
    graph.v_to_e[0].push_back(0);
    graph.v_to_e[0].push_back(5);
    graph.v_to_e[0].push_back(6);
    graph.v_to_e[0].push_back(7);
    graph.v_to_e[0].push_back(18);
    graph.v_to_e[0].push_back(19);
    graph.v_to_e[1].push_back(0);
    graph.v_to_e[1].push_back(1);
    graph.v_to_e[1].push_back(8);
    graph.v_to_e[1].push_back(9);
    graph.v_to_e[2].push_back(1);
    graph.v_to_e[2].push_back(2);
    graph.v_to_e[2].push_back(10);
    graph.v_to_e[2].push_back(11);
    graph.v_to_e[3].push_back(2);
    graph.v_to_e[3].push_back(3);
    graph.v_to_e[3].push_back(6);
    graph.v_to_e[3].push_back(12);
    graph.v_to_e[3].push_back(13);
    graph.v_to_e[3].push_back(20);
    graph.v_to_e[4].push_back(3);
    graph.v_to_e[4].push_back(4);
    graph.v_to_e[4].push_back(14);
    graph.v_to_e[4].push_back(15);
    graph.v_to_e[5].push_back(4);
    graph.v_to_e[5].push_back(5);
    graph.v_to_e[5].push_back(16);
    graph.v_to_e[5].push_back(17);
    graph.v_to_e[6].push_back(7);
    graph.v_to_e[6].push_back(8);
    graph.v_to_e[7].push_back(9);
    graph.v_to_e[7].push_back(10);
    graph.v_to_e[8].push_back(11);
    graph.v_to_e[8].push_back(12);
    graph.v_to_e[9].push_back(13);
    graph.v_to_e[9].push_back(14);
    graph.v_to_e[10].push_back(15);
    graph.v_to_e[10].push_back(16);
    graph.v_to_e[11].push_back(17);
    graph.v_to_e[11].push_back(18);
    graph.v_to_e[12].push_back(19);
    graph.v_to_e[12].push_back(20);

    std::cout << "Testing subgraph path" << std::endl;
    Subgraph sg(21, true);
    sg.edges[0] = false;
    sg.edges[2] = false;
    sg.edges[6] = false;
    sg.edges[8] = false;
    Subgraph path;
    Subgraph path_no_clear(21, false);
    get_subgraph_path(graph, sg, 0, 1, path, path_no_clear);
    for (unsigned i = 0; i < path.edges.size(); ++i)
        std::cout << path.edges[i] << " ";
    std::cout << std::endl;
    for (unsigned i = 0; i < path.vertex_list.size(); ++i)
        std::cout << path.vertex_list[i] << " ";
    std::cout << std::endl << std::endl;

    std::cout << "Testing cycle basis" << std::endl;
    std::deque<Subgraph> basis;
    std::deque<int> pivots;
    std::vector<int> non_pivots;
    unsigned n_fixed = get_cycle_basis(graph, basis, pivots, non_pivots);
    std::cout << "Basis: " << std::endl;
    for (unsigned i = 0; i < basis.size(); ++i) {
        for (unsigned j = 0; j < basis[i].edges.size(); ++j)
            std::cout << basis[i].edges[j] << " ";
        std::cout << std::endl;
        for (unsigned j = 0; j < basis[i].vertex_list.size(); ++j)
            std::cout << basis[i].vertex_list[j] << " ";
        std::cout << std::endl;
    }
    std::cout << "Pivots: " << std::endl;
    for (unsigned i = 0; i < pivots.size(); ++i)
        std::cout << pivots[i] << " ";
    std::cout << std::endl;
    std::cout << "Non-pivots: " << std::endl;
    for (unsigned i = 0; i < non_pivots.size(); ++i)
        std::cout << non_pivots[i] << " ";
    std::cout << std::endl;
    std::cout << "Fixed: " << n_fixed << std::endl << std::endl;

    std::cout << "Testing odd path" << std::endl;
    sg = Subgraph(21, true);
    sg.edges[0] = false;
    sg.edges[5] = false;
    sg.edges[6] = false;
    std::vector<Subgraph> paths;
    get_odd_path(graph, sg, 0, 0, false, -1, paths);
    for (unsigned i = 0; i < paths[0].edges.size(); ++i)
        std::cout << paths[0].edges[i] << " ";
    std::cout << std::endl;
    for (unsigned i = 0; i < paths[0].vertex_list.size(); ++i)
        std::cout << paths[0].vertex_list[i] << " ";
    std::cout << std::endl << std::endl;

    std::cout << "Testing minimize basis" << std::endl;
    std::vector<Subgraph> min_basis;
    minimize_cycle_basis(graph, basis, pivots, n_fixed, non_pivots, min_basis);
    std::cout << "Min basis: " << std::endl;
    for (unsigned i = 0; i < min_basis.size(); ++i) {
        for (unsigned j = 0; j < min_basis[i].edges.size(); ++j)
            std::cout << min_basis[i].edges[j] << " ";
        std::cout << std::endl;
        for (unsigned j = 0; j < min_basis[i].vertex_list.size(); ++j)
            std::cout << min_basis[i].vertex_list[j] << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Testing relevant cycles" << std::endl;
    std::set<Subgraph> relevant_cycles;
    get_relevant_cycles(graph, min_basis, pivots, non_pivots, relevant_cycles);
    std::cout << "Relevant cycles: " << std::endl;
    for (std::set<Subgraph>::iterator iter = relevant_cycles.begin();
           iter != relevant_cycles.end(); ++iter) {
        for (unsigned j = 0; j < iter->edges.size(); ++j)
            std::cout << iter->edges[j] << " ";
        std::cout << std::endl;
        for (unsigned j = 0; j < iter->vertex_list.size(); ++j)
            std::cout << iter->vertex_list[j] << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

static void test_sssr(const char* input_dms) {
    std::cout << "Testing SSSR" << std::endl;
    SystemPtr sys = ImportDMS(input_dms, true);
    MultiIdList rings = GetSSSR(sys, sys->atoms());
    for (unsigned i = 0; i < rings.size(); ++i) {
        for (unsigned j = 0; j < rings[i].size(); ++j)
            std::cout << rings[i][j] << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Testing all relevant cycles" << std::endl;
    rings = GetSSSR(sys, sys->atoms(), true);
    for (unsigned i = 0; i < rings.size(); ++i) {
        for (unsigned j = 0; j < rings[i].size(); ++j)
            std::cout << rings[i][j] << " ";
        std::cout << std::endl;
    }
}

int main(int argc,char **argv){
    test_biconnected_components();
    test_min_cycle_basis();
    for (int i=1; i<argc; i++) {
        test_sssr(argv[i]);
    }
    return 0;
}
