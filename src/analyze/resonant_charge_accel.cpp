#include <boost/numeric/conversion/cast.hpp>
#include <numeric>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#include <bliss/graph.hh>

#include "bond_orders.hxx"

using namespace desres::msys;

namespace lpsolve {
/* for solving the integer linear equations. We put this in a new
   namespace since it doesn't already use one */
#if defined __has_include
#if __has_include(<lp_solve/lp_lib.h>)
// on some systems, the lp_solve headers are dumped right in /usr/lib/include
// or equivalent, rather than all being inside a directory called lp_solve
#include <lp_solve/lp_lib.h>
#else
#include <lp_lib.h>
#endif
#else
// if the compiler doesn't have __has_include, then just hope that the lpsolve
// header is in its own directory.
#include <lp_solve/lp_lib.h>
#endif
} // namespace lpsolve

namespace {
// Given a vector of length N, return a vector of integers of
// length N giving that labels each unique element in the input
// vector with a different integer.
// Example: labelUnique({2.5, 2.5, 3.0, 2.5}) == {0, 0, 1, 0};
template <typename T>
std::vector<int64_t> labelUnique(const std::vector<T> &v_) {
  std::vector<T> v = v_; // make a copy

  // sort and erase duplicates
  std::sort(v.begin(), v.end());
  auto it = std::unique(v.begin(), v.end());
  v.resize(std::distance(v.begin(), it));

  std::vector<int64_t> answer(v_.size());
  for (size_t i = 0; i < v_.size(); i++) {
    auto lower = std::lower_bound(v.begin(), v.end(), v_[i]);
    assert(lower != v.end());
    answer[i] = std::distance(v.begin(), lower);
  }
  return answer;
}

uint32_t safe_shl(uint32_t x, uint8_t y) {
  if (x == 0) {
    return 0;
  }
  assert(y < 32);
  uint32_t z = (1u * x) << y;
  assert((z >> y) == x);
  return z;
}

void blissRecordAutomorphismCallback(void *generators_nvars, unsigned int n,
                                     const unsigned int *aut) {
  static_assert(sizeof(unsigned int) == sizeof(Id));
  std::pair<MultiIdList, uint32_t>* generators_and_nvars = static_cast<std::pair<MultiIdList, uint32_t> *>(generators_nvars);
  assert(generators_nvars != nullptr);
  MultiIdList& gen = generators_and_nvars->first;
  uint32_t nvars = generators_and_nvars->second;
  assert(nvars <= n);

  IdList cycle(nvars);
  for (auto i = 0u; i < nvars; i++) {
    cycle[i] = aut[i];
  }
  gen.push_back(cycle);
}

/*
 * Generate all permuations of vector `base` under the permutation group
 * spanned by `generators`.
 *
 * For example, if generators is {{1, 2, 3, 0}, {0, 2, 1, 3}}, these are
 * two permutations written in array notation. The first takes elements
 * 0->1, 1->2, 2->3, and 3->0, a cyclic permutation. The second swaps
 * elements 2 and 3. Arbitrary compositions and powers of these generators
 * give rise to a whole family of possible permutations. In fact, the
 * compositions of these 2 permutations happen to span all possible
 * permutations on 4 elements.
 *
 * Applying these permutations to the vector {"a", "b", "c", "d"} yields
 * all 24 permutations. Applying them to the vector {"a", "a", "b", "c"}
 * would yield 12. Applying only the permutation {1, 2, 3, 0} to the vector
 * {"a", "b", "c", "d"} would yield the 4 cyclic permutations. And so on.
 *
 * As some more example, applying the generators {{1, 2, 3, 0}, {1, 0, 3,
 * 2}} to {"a", "b", "c", "d"} gives 8 permutations, and applying those
 * generators to {"a", "a", "b", "b"} would give 4 permutations.
 *
 */
template <typename T>
std::vector<std::vector<T>>
make_permutations(const std::vector<T> &base,
                  const std::vector<std::vector<uint32_t>> &generators) {

  for (const auto &g : generators) {
    uint64_t sum = 0;
    for (const auto &p : g) {
      assert(p >= 0 && p < g.size());
      sum += p;
    }
    assert(sum == g.size() * (g.size() - 1) / 2);
    assert(g.size() == base.size());
  }

  auto hash = [](const std::vector<T> &v) {
    std::hash<T> hasher;
    size_t seed = 0;
    for (auto const &i : v) {
      seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  };

  auto apply = [](const std::vector<T> &a,
                  const std::vector<uint32_t> &gen) -> std::vector<T> {
    std::vector<T> ag(gen.size());
    for (size_t i = 0u; i < gen.size(); i++) {
      ag[i] = a[gen[i]];
    }
    return ag;
  };

  std::unordered_set<std::vector<T>, decltype(hash)> expanded{{}, hash};
  std::vector<std::vector<T>> current{{base}};
  std::vector<std::vector<T>> next;

  while (current.size() > 0) {
    for (const auto &elem : current) {
      for (const auto &g : generators) {
        auto ag = apply(elem, g);
        if (expanded.find(ag) == expanded.end()) {
          next.push_back(ag);
        }
      }
      expanded.insert(elem);
    }
    std::swap(current, next);
    next.clear();
  }

  return std::vector<std::vector<T>>{expanded.begin(), expanded.end()};
}

} // namespace

namespace desres {
namespace msys {

std::vector<std::vector<int>> ComponentAssigner::generate_resonance_ilp_permutations(
    lpsolve::_lprec* lp,
    const std::vector<int>& solution
) {
  bliss::Graph lpGraph(0);

  auto makeBlissColor = [&](int8_t category, int64_t color) {
    assert(category < 4);
    assert(color >= 0 && color < UINT32_MAX);
    return category + safe_shl(static_cast<uint32_t>(color), 2);
  };

  std::vector<int> orig2presolved(1 + lpsolve::get_Norig_columns(lp), -1);
  for (int i = 1; i < 1 + lpsolve::get_Ncolumns(lp); i++) {
    auto iorig = lpsolve::get_orig_index(lp, lpsolve::get_Nrows(lp) + i);
    assert (iorig >= 0 && iorig < 1 + lpsolve::get_Norig_columns(lp));
    orig2presolved[iorig] = i;
  }

  // Add one node per variable in the ILP and color them
  // by unique (objective coefficient and upper/lower bound)
  // in the ILP objective function
  // (Each of these atoms get atomic_number 1 to distinguish
  // them from nodes in other categories)
  {
    std::vector<std::tuple<double, double, double, size_t>> coefVector;
    for (int jj = 0; jj <= lpsolve::get_Norig_columns(lp); jj++) {
      auto j = orig2presolved[jj];
      if (j >= 0) {
        coefVector.push_back({
            (j > 0) ? lpsolve::get_mat(lp, 0, j) : 0,
            (j > 0) ? lpsolve::get_lowbo(lp, j) : 0,
            (j > 0) ? lpsolve::get_upbo(lp, j) : 0,
            (j > 0) ? lpsolve::is_int(lp, j) : 0,
        });
      } else {
        // presolved, so give it something unique
        auto u = std::hash<int>()(jj);
        coefVector.push_back({0, 0, 0, u});
      }
    }

    std::vector<int64_t> coefColors = labelUnique(coefVector);
    assert(coefColors.size() == coefVector.size());

    for (size_t j = 0; j < coefColors.size(); j++) {
      lpGraph.add_vertex(makeBlissColor(1, coefColors[j]));
    }
  }

  // Add one node per constraint in the ILP and color them
  // by (RHS coefficient, the type of the constraint, and any
  // "range", which is when the row has both an upper and lower bound)
  // (Each of these atoms get atomic_number 2 to distinguish
  // them from nodes in other categories)
  {
    std::vector<std::tuple<double, int, double>> rhsVector;
    for (auto i = 1; i <= lpsolve::get_Nrows(lp); i++) {
      rhsVector.push_back({lpsolve::get_rh(lp, i),
                           lpsolve::get_constr_type(lp, i),
                           lpsolve::get_rh_range(lp, i)});
    }
    std::vector<int64_t> rhsColors = labelUnique(rhsVector);
    assert(rhsColors.size() == rhsVector.size());
    for (size_t i = 0; i < rhsColors.size(); i++) {
      lpGraph.add_vertex(makeBlissColor(2, rhsColors[i]));
    }
  }

  // Add one node per nonzero element in the constraint matrix and
  // color them by the matrix element
  {
    std::vector<double> linearizedConstraintMatrixElements;
    std::vector<double> row(1 + lpsolve::get_Ncolumns(lp));
    for (auto i = 1; i <= lpsolve::get_Nrows(lp); i++) {
      lpsolve::get_row(lp, i, &row[0]);
      for (auto el : row) {
        linearizedConstraintMatrixElements.push_back(el);
      }
    }
    std::vector<int64_t> constraintMatrixColors =
        labelUnique(linearizedConstraintMatrixElements);

    for (auto i = 1; i <= lpsolve::get_Nrows(lp); i++) {
      lpsolve::get_row(lp, i, &row[0]);
      for (size_t j = 0; j < row.size(); j++) {
        if (row[j] == 0) {
          continue;
        }

        // get the ids associated with the (i, j) position in the matrix.
        // there are three node ids:
        // (a) the node associated with the variable which is just j, except
        //     that we need to convert from post- pre-solved numbering to
        //     original numbering.
        // (b) the node associated with the constraint. since these were
        // inserted after the variable nodes, theyy're offset
        // (c) a new node for the connecting matrix element

        int variableNodeId = lpsolve::get_orig_index(lp, lpsolve::get_Nrows(lp) + j);
        assert(variableNodeId >= 0 && variableNodeId < 1 + lpsolve::get_Norig_columns(lp));
        msys::Id constraintNodeId = lpsolve::get_Ncolumns(lp) + i;
        // assert (lpGraph->atom(variableNodeId).atomic_number == 1);
        // assert (lpGraph->atom(constraintNodeId).atomic_number == 2);

        auto linearMatrixPosition = ((i - 1) * row.size()) + j;
        assert(row[j] ==
               linearizedConstraintMatrixElements[linearMatrixPosition]);

        auto connectorNodeId = lpGraph.add_vertex(
            makeBlissColor(3, constraintMatrixColors[linearMatrixPosition]));

        lpGraph.add_edge(variableNodeId, connectorNodeId);
        lpGraph.add_edge(connectorNodeId, constraintNodeId);
      }
    }
  }

  bliss::Stats stats;
  std::pair<MultiIdList, uint32_t> generators_and_vars;
  generators_and_vars.second = 1 + lpsolve::get_Norig_columns(lp);
  lpGraph.find_automorphisms(stats, &blissRecordAutomorphismCallback,
                             &generators_and_vars);
  const auto& generators = generators_and_vars.first;

  auto solution_permutations = make_permutations<int>(
      solution,
      generators
  );

  return solution_permutations;
}

} // namespace msys
} // namespace desres
