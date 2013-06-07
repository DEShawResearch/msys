#include "atommatch.hxx"
#include "io.hxx"

#include <iostream>

using namespace desres::msys;
using namespace desres::fep_atommatch;
using namespace desres::fep_atommatch::helpers;

void test_init_graph() {
    SystemPtr mol = Load("/proj/desres/root/Linux/x86_64/dms_inputs/1.5.4/share/ww.dms", true);
    MultiIdList frags;
    unsigned num = mol->updateFragids(&frags);
    for (unsigned i = 0; i < num; ++i) {
        if (i == 1) continue;
        BOOST_FOREACH(Id atom, frags[i])
            mol->delAtom(atom);
    }
    std::vector<int> atom_idx;
    GraphRepr g;
    init_graph(mol, atom_idx, g);
    std::cout << "Atom idx: ";
    for (unsigned i = 0; i < atom_idx.size(); ++i)
        std::cout << atom_idx[i] << " ";
    std::cout << std::endl;
    std::cout << "Graph: " << std::endl;
    for (unsigned i = 0; i < g.v_to_e.size(); ++i) {
        std::cout << "Vertex " << i << ": ";
        for (unsigned j = 0; j < g.v_to_e[i].size(); ++j)
            std::cout << g.v_to_e[i][j] << " ";
        std::cout << std::endl;
    }
    for (unsigned i = 0; i < g.edges.size(); ++i)
        std::cout << "Edge " << i << ": " << g.edges[i].first << " " << g.edges[i].second << std::endl;
}

void test_partition_graph() {
    GraphRepr g;
    g.v_to_e.resize(20);
    g.edges.resize(24);
    g.v_to_e[0].push_back(0);
    g.v_to_e[0].push_back(5);
    g.v_to_e[1].push_back(0);
    g.v_to_e[1].push_back(1);
    g.v_to_e[2].push_back(1);
    g.v_to_e[2].push_back(2);
    g.v_to_e[3].push_back(2);
    g.v_to_e[3].push_back(3);
    g.v_to_e[4].push_back(3);
    g.v_to_e[4].push_back(4);
    g.v_to_e[5].push_back(4);
    g.v_to_e[5].push_back(5);
    g.v_to_e[5].push_back(6);
    g.v_to_e[5].push_back(9);
    g.v_to_e[6].push_back(6);
    g.v_to_e[6].push_back(7);
    g.v_to_e[6].push_back(21);
    g.v_to_e[7].push_back(7);
    g.v_to_e[7].push_back(8);
    g.v_to_e[7].push_back(10);
    g.v_to_e[7].push_back(22);
    g.v_to_e[8].push_back(8);
    g.v_to_e[8].push_back(9);
    g.v_to_e[8].push_back(11);
    g.v_to_e[9].push_back(10);
    g.v_to_e[9].push_back(11);
    g.v_to_e[9].push_back(12);
    g.v_to_e[10].push_back(12);
    g.v_to_e[10].push_back(13);
    g.v_to_e[10].push_back(23);
    g.v_to_e[11].push_back(13);
    g.v_to_e[11].push_back(14);
    g.v_to_e[12].push_back(14);
    g.v_to_e[12].push_back(15);
    g.v_to_e[12].push_back(23);
    g.v_to_e[13].push_back(15);
    g.v_to_e[13].push_back(16);
    g.v_to_e[14].push_back(16);
    g.v_to_e[14].push_back(17);
    g.v_to_e[15].push_back(17);
    g.v_to_e[15].push_back(18);
    g.v_to_e[15].push_back(20);
    g.v_to_e[16].push_back(18);
    g.v_to_e[16].push_back(19);
    g.v_to_e[17].push_back(19);
    g.v_to_e[17].push_back(20);
    g.v_to_e[18].push_back(21);
    g.v_to_e[19].push_back(22);
    g.edges[0] = Edge(0,1);
    g.edges[1] = Edge(1,2);
    g.edges[2] = Edge(2,3);
    g.edges[3] = Edge(3,4);
    g.edges[4] = Edge(4,5);
    g.edges[5] = Edge(0,5);
    g.edges[6] = Edge(5,6);
    g.edges[7] = Edge(6,7);
    g.edges[8] = Edge(7,8);
    g.edges[9] = Edge(5,8);
    g.edges[10] = Edge(7,9);
    g.edges[11] = Edge(8,9);
    g.edges[12] = Edge(9,10);
    g.edges[13] = Edge(10,11);
    g.edges[14] = Edge(11,12);
    g.edges[15] = Edge(12,13);
    g.edges[16] = Edge(13,14);
    g.edges[17] = Edge(14,15);
    g.edges[18] = Edge(15,16);
    g.edges[19] = Edge(16,17);
    g.edges[20] = Edge(15,17);
    g.edges[21] = Edge(6,18);
    g.edges[22] = Edge(7,19);
    g.edges[23] = Edge(10,12);
    std::vector<GraphRepr> components;
    std::vector<std::vector<int> > components_idx;
    GraphRepr tree;
    TreeEdges tree_edges_biconnected;
    partition_graph_biconnected(g, components, components_idx, tree, tree_edges_biconnected);
    for (unsigned i = 0; i < components_idx.size(); ++i) {
        std::cout << "Component " << i << ": " << std::endl;
        for (unsigned j = 0; j < components[i].v_to_e.size(); ++j) {
            std::cout << "Vertex " << j << "(" << components_idx[i][j] << "): ";
            BOOST_FOREACH(int edge, components[i].v_to_e[j])
                std::cout << edge << " ";
            std::cout << std::endl;
        }
        for (unsigned j = 0; j < components[i].edges.size(); ++j)
            std::cout << "Edge " << j << ": " << components[i].edges[j].first << " " << components[i].edges[j].second << std::endl;
    }
    std::cout << "Tree: " << std::endl;
    for (unsigned i = 0; i < tree.v_to_e.size(); ++i) {
        std::cout << "Vertex " << i << ": ";
        for (unsigned j = 0; j < tree.v_to_e[i].size(); ++j)
            std::cout << tree.v_to_e[i][j] << " ";
        std::cout << std::endl;
    }
    for (unsigned i = 0; i < tree.edges.size(); ++i)
        std::cout << "Edge " << i << ": " << tree.edges[i].first << " " << tree.edges[i].second << std::endl;
    std::cout << "Tree edges: " << std::endl;
    for (unsigned i = 0; i < tree_edges_biconnected.size(); ++i) {
        for (unsigned j = 0; j < tree_edges_biconnected[i].size(); ++j) {
            for (unsigned k = 0; k < tree_edges_biconnected[i][j].size(); ++k) {
                std::cout << "(" << i << "," << j << ")~(" << tree_edges_biconnected[i][j][k].first << "," << tree_edges_biconnected[i][j][k].second << ") " << std::endl;
            }
        }
    }
}

double rep1(SystemPtr mol1, const IdList& atoms1,
        SystemPtr mol2, const IdList& atoms2) {
    return 1;
}

double rep2(SystemPtr mol1, const IdList& atoms1,
        SystemPtr mol2, const IdList& atoms2) {
    if (atoms1.size() == 2)
        return 1;
    if (atoms1.size() != atoms2.size())
        FAIL("BUG---Atom lists of different sizes");
    for (unsigned i = 0; i < atoms1.size(); ++i) {
        if (mol1->bondsForAtom(atoms1[i]).size()
                != mol2->bondsForAtom(atoms2[i]).size())
            return -1;
    }
    return 1;
}

double rep3(SystemPtr mol1, const IdList& atoms1,
        SystemPtr mol2, const IdList& atoms2) {
    if (atoms1.size() != atoms2.size())
        FAIL("BUG---Atom lists of different sizes");
    for (unsigned i = 0; i < atoms1.size(); ++i) {
        if (mol1->bondsForAtom(atoms1[i]).size()
                != mol2->bondsForAtom(atoms2[i]).size())
            return -1;
    }
    return 1;
}

double rep4(SystemPtr mol1, const IdList& atoms1,
        SystemPtr mol2, const IdList& atoms2) {
    return 0;
}

void test_isomorphisms() {
    SystemPtr mol1 = System::create();
    Id res1 = mol1->addResidue(mol1->addChain());
    SystemPtr mol2 = System::create();
    Id res2 = mol2->addResidue(mol2->addChain());
    for (unsigned i = 0; i < 6; ++i) {
        mol1->addAtom(res1);
        mol2->addAtom(res2);
    }
    for (unsigned i = 0; i < 6; ++i) {
        mol1->addBond(i, (i+1)%6);
        mol2->addBond(i, (i+1)%6);
    }
    for (unsigned i = 0; i < 6; ++i) {
        Id a = mol1->addAtom(res1);
        mol1->addBond(i, a);
        Id b = mol2->addAtom(res2);
        Id c = mol2->addAtom(res2);
        mol2->addBond(i, b);
        mol2->addBond(i, c);
    }
    std::vector<int> atom_idx1;
    std::vector<int> atom_idx2;
    GraphRepr g1;
    GraphRepr g2;
    init_graph(mol1, atom_idx1, g1);
    init_graph(mol2, atom_idx2, g2);
    std::vector<GraphRepr> components1;
    std::vector<GraphRepr> components2;
    std::vector<std::vector<int> > components_idx1;
    std::vector<std::vector<int> > components_idx2;
    GraphRepr tree1;
    GraphRepr tree2;
    /*
    std::vector<std::vector<std::pair<int, int> > > tree_edges1;
    std::vector<std::vector<std::pair<int, int> > > tree_edges2;
    partition_graph_biedgeconnected(g1, components1, components_idx1, tree1, tree_edges1);
    partition_graph_biedgeconnected(g2, components2, components_idx2, tree2, tree_edges2);
    */
    TreeEdges tree_edges1;
    TreeEdges tree_edges2;
    partition_graph_biconnected(g1, components1, components_idx1, tree1, tree_edges1);
    partition_graph_biconnected(g2, components2, components_idx2, tree2, tree_edges2);
    unsigned ind1 = 0;
    for (; ind1 < components1.size(); ++ind1)
        if (components1[ind1].v_to_e.size() == 6) break;
    unsigned ind2 = 0;
    for (; ind2 < components2.size(); ++ind2)
        if (components2[ind2].v_to_e.size() == 6) break;
    std::cout << "Graph 1: " << std::endl;
    for (unsigned i = 0; i < components1[ind1].v_to_e.size(); ++i) {
        std::cout << "Vertex " << i << ": ";
        for (unsigned j = 0; j < components1[ind1].v_to_e[i].size(); ++j)
            std::cout << components1[ind1].v_to_e[i][j] << " ";
        std::cout << std::endl;
    }
    for (unsigned i = 0; i < components1[ind1].edges.size(); ++i)
        std::cout << "Edge " << i << ": " << components1[ind1].edges[i].first << " " << components1[ind1].edges[i].second << std::endl;
    std::cout << "Graph 2: " << std::endl;
    for (unsigned i = 0; i < components2[ind2].v_to_e.size(); ++i) {
        std::cout << "Vertex " << i << ": ";
        for (unsigned j = 0; j < components2[ind2].v_to_e[i].size(); ++j)
            std::cout << components2[ind2].v_to_e[i][j] << " ";
        std::cout << std::endl;
    }
    for (unsigned i = 0; i < components2[ind2].edges.size(); ++i)
        std::cout << "Edge " << i << ": " << components2[ind2].edges[i].first << " " << components2[ind2].edges[i].second << std::endl;
    Isomorphisms isos1;
    Isomorphisms isos2;
    isomorphisms(components1[ind1], components_idx1[ind1], atom_idx1, mol1,
            components2[ind2], components_idx2[ind2], atom_idx2, mol2,
            ScoreFctPtr(new ScoreFctC(rep1)), isos1);
    isomorphisms(components1[ind1], components_idx1[ind1], atom_idx1, mol1,
            components2[ind2], components_idx2[ind2], atom_idx2, mol2,
            ScoreFctPtr(new ScoreFctC(rep2)), isos2);
    std::cout << "Isomorphisms 1: " << std::endl;
    for (unsigned i = 0; i < isos1.size(); ++i) {
        for (unsigned j = 0; j < isos1[i].perm.size(); ++j)
            std::cout << isos1[i].perm[j] << " ";
        std::cout << "; score " << isos1[i].score << std::endl;
    }
    std::cout << "Isomorphisms 2: " << std::endl;
    for (unsigned i = 0; i < isos2.size(); ++i) {
        for (unsigned j = 0; j < isos2[i].perm.size(); ++j)
            std::cout << isos2[i].perm[j] << " ";
        std::cout << "; score " << isos2[i].score << std::endl;
    }
}

void test_get_levels() {
    GraphRepr g;
    g.v_to_e.resize(8);
    g.edges.resize(7);
    g.v_to_e[0].push_back(0);
    g.v_to_e[0].push_back(1);
    g.v_to_e[1].push_back(0);
    g.v_to_e[1].push_back(2);
    g.v_to_e[1].push_back(3);
    g.v_to_e[2].push_back(1);
    g.v_to_e[3].push_back(2);
    g.v_to_e[4].push_back(3);
    g.v_to_e[4].push_back(4);
    g.v_to_e[4].push_back(6);
    g.v_to_e[5].push_back(4);
    g.v_to_e[5].push_back(5);
    g.v_to_e[6].push_back(5);
    g.v_to_e[7].push_back(6);
    g.edges[0] = Edge(0,1);
    g.edges[1] = Edge(0,2);
    g.edges[2] = Edge(1,3);
    g.edges[3] = Edge(1,4);
    g.edges[4] = Edge(4,5);
    g.edges[5] = Edge(5,6);
    g.edges[6] = Edge(4,7);
    std::vector<std::vector<int> > levels1;
    std::vector<int> level_map1;
    std::vector<std::vector<int> > levels2;
    std::vector<int> level_map2;
    std::vector<std::vector<int> > levels3;
    std::vector<int> level_map3;
    get_levels(g, 0, levels1, level_map1, false);
    get_levels(g, 4, levels2, level_map2, false);
    get_levels(g, 4, levels3, level_map3, true);
    std::cout << "Run 1: " << std::endl;
    std::cout << "Level map: ";
    BOOST_FOREACH(int level, level_map1)
        std::cout << level << " ";
    std::cout << std::endl;
    std::cout << "Levels: ";
    BOOST_FOREACH(const std::vector<int>& levels, levels1) {
        std::cout << "(";
        BOOST_FOREACH(int l, levels)
            std::cout << l << " ";
        std::cout << ") ";
    }
    std::cout << std::endl;
    std::cout << "Run 2: " << std::endl;
    std::cout << "Level map: ";
    BOOST_FOREACH(int level, level_map2)
        std::cout << level << " ";
    std::cout << std::endl;
    std::cout << "Levels: ";
    BOOST_FOREACH(const std::vector<int>& levels, levels2) {
        std::cout << "(";
        BOOST_FOREACH(int l, levels)
            std::cout << l << " ";
        std::cout << ") ";
    }
    std::cout << std::endl;
    std::cout << "Run 3: " << std::endl;
    std::cout << "Level map: ";
    BOOST_FOREACH(int level, level_map3)
        std::cout << level << " ";
    std::cout << std::endl;
    std::cout << "Levels: ";
    BOOST_FOREACH(const std::vector<int>& levels, levels3) {
        std::cout << "(";
        BOOST_FOREACH(int l, levels)
            std::cout << l << " ";
        std::cout << ") ";
    }
    std::cout << std::endl;
}

void test_get_perms() {
    std::vector<std::vector<int> > perms1;
    std::vector<std::vector<int> > perms2;
    std::vector<std::vector<int> > perms3;
    std::vector<std::vector<int> > perms4;
    get_perms(2, 4, perms1);
    get_perms(2, 4, perms2);
    get_perms(4, 2, perms3);
    get_perms(3, 3, perms4);
    std::cout << "Perms 1: " << std::endl;
    BOOST_FOREACH(const std::vector<int>& perm, perms1) {
        BOOST_FOREACH(int p, perm)
            std::cout << p << " ";
        std::cout << std::endl;
    }
    std::cout << "Perms 2: " << std::endl;
    BOOST_FOREACH(const std::vector<int>& perm, perms2) {
        BOOST_FOREACH(int p, perm)
            std::cout << p << " ";
        std::cout << std::endl;
    }
    std::cout << "Perms 3: " << std::endl;
    BOOST_FOREACH(const std::vector<int>& perm, perms3) {
        BOOST_FOREACH(int p, perm)
            std::cout << p << " ";
        std::cout << std::endl;
    }
    std::cout << "Perms 4: " << std::endl;
    BOOST_FOREACH(const std::vector<int>& perm, perms4) {
        BOOST_FOREACH(int p, perm)
            std::cout << p << " ";
        std::cout << std::endl;
    }
}

void test_options_append() {
    std::vector<DPNodes> options(1, DPNodes());
    std::vector<int> inds;
    inds.push_back(0);
    inds.push_back(1);
    std::cout << "Iter 0: " << std::endl;
    BOOST_FOREACH(const DPNodes& nodes, options) {
        BOOST_FOREACH(const DPNode& node, nodes)
            std::cout << "(" << node.v1 << "," << node.v2 << "," << node.iso_ind << ") ";
        std::cout << std::endl;
    }
    options_append(options, 0, 0, inds);
    std::cout << "Iter 1: " << std::endl;
    BOOST_FOREACH(const DPNodes& nodes, options) {
        BOOST_FOREACH(const DPNode& node, nodes)
            std::cout << "(" << node.v1 << "," << node.v2 << "," << node.iso_ind << ") ";
        std::cout << std::endl;
    }
    options_append(options, 1, 1, inds);
    std::cout << "Iter 2: " << std::endl;
    BOOST_FOREACH(const DPNodes& nodes, options) {
        BOOST_FOREACH(const DPNode& node, nodes)
            std::cout << "(" << node.v1 << "," << node.v2 << "," << node.iso_ind << ") ";
        std::cout << std::endl;
    }
    options_append(options, 2, 2, inds);
    std::cout << "Iter 3: " << std::endl;
    BOOST_FOREACH(const DPNodes& nodes, options) {
        BOOST_FOREACH(const DPNode& node, nodes)
            std::cout << "(" << node.v1 << "," << node.v2 << "," << node.iso_ind << ") ";
        std::cout << std::endl;
    }
}

void test_dp_update(int root1, int root2, ScoreFctPtr rep) {
    SystemPtr mol1 = System::create();
    Id res1 = mol1->addResidue(mol1->addChain());
    SystemPtr mol2 = System::create();
    Id res2 = mol2->addResidue(mol2->addChain());
    for (unsigned i = 0; i < 6; ++i) {
        mol1->addAtom(res1);
        mol2->addAtom(res2);
    }
    for (unsigned i = 0; i < 6; ++i) {
        mol1->addBond(i, (i+1)%6);
        mol2->addBond(i, (i+1)%6);
    }
    for (unsigned i = 0; i < 6; ++i) {
        Id a = mol1->addAtom(res1);
        mol1->addBond(i, a);
        Id b = mol2->addAtom(res2);
        Id c = mol2->addAtom(res2);
        mol2->addBond(i, b);
        mol2->addBond(i, c);
    }
    std::vector<int> atom_idx1;
    std::vector<int> atom_idx2;
    GraphRepr g1;
    GraphRepr g2;
    init_graph(mol1, atom_idx1, g1);
    init_graph(mol2, atom_idx2, g2);
    std::vector<GraphRepr> components1;
    std::vector<GraphRepr> components2;
    std::vector<std::vector<int> > components_idx1;
    std::vector<std::vector<int> > components_idx2;
    GraphRepr tree1;
    GraphRepr tree2;
    /*
    std::vector<std::vector<std::pair<int, int> > > tree_edges1;
    std::vector<std::vector<std::pair<int, int> > > tree_edges2;
    partition_graph_biedgeconnected(g1, components1, components_idx1, tree1, tree_edges1);
    partition_graph_biedgeconnected(g2, components2, components_idx2, tree2, tree_edges2);
    */
    TreeEdges tree_edges1;
    TreeEdges tree_edges2;
    partition_graph_biconnected(g1, components1, components_idx1, tree1, tree_edges1);
    partition_graph_biconnected(g2, components2, components_idx2, tree2, tree_edges2);
    std::vector<std::vector<int> > levels1;
    std::vector<int> level_map1;
    std::vector<std::vector<int> > levels2;
    std::vector<int> level_map2;
    get_levels(tree1, root1, levels1, level_map1, false);
    get_levels(tree2, root2, levels2, level_map2, false);
    std::cout << "Levels 1: ";
    BOOST_FOREACH(const std::vector<int>& levels, levels1) {
        std::cout << "(";
        BOOST_FOREACH(int l, levels)
            std::cout << l << "[" << components_idx1[l].size() << "] ";
        std::cout << ") ";
    }
    std::cout << std::endl;
    std::cout << "Levels 2: ";
    BOOST_FOREACH(const std::vector<int>& levels, levels2) {
        std::cout << "(";
        BOOST_FOREACH(int l, levels)
            std::cout << l << "[" << components_idx2[l].size() << "] ";
        std::cout << ") ";
    }
    std::cout << std::endl;
    DPMat best_matches(tree1.v_to_e.size(), std::vector<std::vector<DPBest> >(
                tree2.v_to_e.size()));
    IsoMat iso_mat(components1.size(),
            std::vector<Isomorphisms>(components2.size()));
    for (unsigned i = 0; i < components1.size(); ++i)
        for (unsigned j = 0; j < components2.size(); ++j)
            isomorphisms(components1[i], components_idx1[i], atom_idx1,
                    mol1, components2[j], components_idx2[j], atom_idx2, mol2,
                    rep, iso_mat[i][j]);
    for (int j = levels1.size() - 1; j >= 0; --j) {
        if (j >= int(levels2.size())) continue;
        BOOST_FOREACH(int v1, levels1[j]) {
            BOOST_FOREACH(int v2, levels2[j]) {
                if (iso_mat[v1][v2].size() == 0)
                    continue;
                dp_update(v1, v2, best_matches, tree_edges1[v1], level_map1,
                        tree_edges2[v2], level_map2, iso_mat);
                for (unsigned i = 0; i < best_matches[v1][v2].size(); ++i) {
                    std::cout << "Best matches " << v1 << " " << v2 << " " << i << ": " << best_matches[v1][v2][i].score << std::endl;
                    BOOST_FOREACH(const DPNodes& option, best_matches[v1][v2][i].options) {
                        BOOST_FOREACH(const DPNode& node, option)
                            std::cout << "(" << node.v1 << "," << node.v2 << "," << node.iso_ind << ") ";
                        std::cout << std::endl;
                    }
                }
            }
        }
    }
}

void test_fep_atommatch(ScoreFctPtr rep) {
    SystemPtr mol1 = System::create();
    Id res1 = mol1->addResidue(mol1->addChain());
    SystemPtr mol2 = System::create();
    Id res2 = mol2->addResidue(mol2->addChain());
    for (unsigned i = 0; i < 4; ++i) {
        mol1->addAtom(res1);
        mol2->addAtom(res2);
    }
    mol1->atom(0).x = 1;
    mol1->atom(0).y = 0;
    mol1->atom(0).z = 4;
    mol1->atom(1).x = 0;
    mol1->atom(1).y = 1;
    mol1->atom(1).z = 4;
    mol1->atom(2).x = -1;
    mol1->atom(2).y = 0;
    mol1->atom(2).z = 4;
    mol1->atom(3).x = 0;
    mol1->atom(3).y = -1;
    mol1->atom(3).z = 4;
    mol2->atom(0).x = 1;
    mol2->atom(0).y = 0;
    mol2->atom(0).z = 0;
    mol2->atom(1).x = 0;
    mol2->atom(1).y = 1;
    mol2->atom(1).z = 0;
    mol2->atom(2).x = -1;
    mol2->atom(2).y = 0;
    mol2->atom(2).z = 0;
    mol2->atom(3).x = 0;
    mol2->atom(3).y = -1;
    mol2->atom(3).z = 0;
    for (unsigned i = 0; i < 4; ++i) {
        mol1->addBond(i, (i+1)%4);
        mol2->addBond(i, (i+1)%4);
    }
    for (unsigned i = 0; i < 4; ++i) {
        Id a = mol1->addAtom(res1);
        mol1->addBond(i, a);
        Id b = mol2->addAtom(res2);
        Id c = mol2->addAtom(res2);
        mol2->addBond(i, b);
        mol2->addBond(i, c);
    }
    mol1->atom(4).x = 2;
    mol1->atom(4).y = 0;
    mol1->atom(4).z = 4;
    mol1->atom(5).x = 0;
    mol1->atom(5).y = 2;
    mol1->atom(5).z = 4;
    mol1->atom(6).x = -2;
    mol1->atom(6).y = 0;
    mol1->atom(6).z = 4;
    mol1->atom(7).x = 0;
    mol1->atom(7).y = -2;
    mol1->atom(7).z = 4;
    mol2->atom(4).x = 1;
    mol2->atom(4).y = 0;
    mol2->atom(4).z = 1;
    mol2->atom(6).x = 0;
    mol2->atom(6).y = 1;
    mol2->atom(6).z = 1;
    mol2->atom(8).x = -1;
    mol2->atom(8).y = 0;
    mol2->atom(8).z = 1;
    mol2->atom(10).x = 0;
    mol2->atom(10).y = -1;
    mol2->atom(10).z = 1;
    mol2->atom(5).x = 2;
    mol2->atom(5).y = 0;
    mol2->atom(5).z = 0;
    mol2->atom(7).x = 0;
    mol2->atom(7).y = 2;
    mol2->atom(7).z = 0;
    mol2->atom(9).x = -2;
    mol2->atom(9).y = 0;
    mol2->atom(9).z = 0;
    mol2->atom(11).x = 0;
    mol2->atom(11).y = -2;
    mol2->atom(11).z = 0;
    MatchList match;
    FEPAtomMatch(mol2, mol1, rep, match);
    for (unsigned i = 0; i < match.size(); ++i)
        std::cout << "(" << match[i].first << "," << match[i].second << ") ";
    std::cout << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "TESTING INIT GRAPH" << std::endl;
    test_init_graph();
    std::cout << "TESTING PARTITION GRAPH" << std::endl;
    test_partition_graph();
    std::cout << "TESTING ISOMORPHISMS" << std::endl;
    test_isomorphisms();
    std::cout << "TESTING GET LEVELS" << std::endl;
    test_get_levels();
    std::cout << "TESTING GET PERMS" << std::endl;
    test_get_perms();
    std::cout << "TESTING OPTIONS APPEND" << std::endl;
    test_options_append();
    std::cout << "TESTING DP UPDATE REP1" << std::endl;
    test_dp_update(1, 3, ScoreFctPtr(new ScoreFctC(rep1)));
    std::cout << "TESTING DP UPDATE REP2" << std::endl;
    test_dp_update(1, 3, ScoreFctPtr(new ScoreFctC(rep2)));
    std::cout << "TESTING DP UPDATE REP3" << std::endl;
    test_dp_update(1, 3, ScoreFctPtr(new ScoreFctC(rep3)));
    std::cout << "TESTING DP UPDATE REP4" << std::endl;
    test_dp_update(1, 4, ScoreFctPtr(new ScoreFctC(rep4)));
    std::cout << "TESTING FEP ATOMMATCH REP1" << std::endl;
    test_fep_atommatch(ScoreFctPtr(new ScoreFctC(rep1)));
    std::cout << "TESTING FEP ATOMMATCH REP2" << std::endl;
    test_fep_atommatch(ScoreFctPtr(new ScoreFctC(rep2)));
    std::cout << "TESTING FEP ATOMMATCH REP3" << std::endl;
    test_fep_atommatch(ScoreFctPtr(new ScoreFctC(rep3)));
    std::cout << "TESTING FEP ATOMMATCH REP4" << std::endl;
    test_fep_atommatch(ScoreFctPtr(new ScoreFctC(rep4)));
    return 0;
}
