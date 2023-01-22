#include "circuits/resistor_network.hpp"

#include "datastructures/disjoint_set_forest.hpp"
#include "equal.hpp"
#include "linmath/contiguous_matrix.hpp"
#include "linmath/linear_solver.hpp"
#include "linmath/matrix.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

#include <iterator>
#include <utility>

#include <unordered_map>
#include <unordered_set>

#include <boost/functional/hash.hpp>

namespace throttle::circuits {

namespace detail {

connected_resistor_network::connected_resistor_network(circuit_graph_type graph) : m_graph{graph} {
  for (const auto &v : m_graph) {
    for (const auto &p : v.second) {
      auto first = v.first, second = p.first;
      auto [res, emf] = p.second;
      if (first < second && throttle::is_roughly_equal(res, 0.0)) m_short_circuits.push_back({first, second, emf});
    }
  }
}

solution connected_resistor_network::solve() const {
  if (m_graph.empty()) throw std::invalid_argument{"Network can't be empty"};

  // Maps indexes 0, 1, .... to corrensponding iterators in the unordered map that represents the input.
  std::unordered_map<unsigned, decltype(m_graph)::const_iterator> iterator_map;
  // Maps indexes from input to 0, 1, ....
  std::unordered_map<unsigned, unsigned> index_map;
  // Maps pairs of input indexes to current variable
  std::unordered_map<std::pair<unsigned, unsigned>, unsigned, boost::hash<std::pair<unsigned, unsigned>>>
      short_circuit_current_map;

  unsigned j = 0;
  for (auto start = std::next(m_graph.begin()), end = m_graph.end(); start != end; ++start, ++j) {
    iterator_map[j] = start;
    index_map[start->first] = j;
  }

  const auto size = iterator_map.size();
  const auto num_short_circuits = m_short_circuits.size();

  for (j = size; const auto &v : m_short_circuits) {
    short_circuit_current_map.insert({{v.first, v.second}, j});
    ++j;
  }

  const auto make_extended_system = [&]() {
    const auto sz = size + num_short_circuits;

    linmath::contiguous_matrix_d extended_matrix{sz, sz + 1};

    for (const auto &v : iterator_map) {
      const auto &[index, map_iter] = v;
      auto row = extended_matrix[index];

      for (const auto &a : map_iter->second) {
        auto [res, emf] = a.second;

        if (throttle::is_roughly_equal(res, 0.0)) {
          const auto initial_index = iterator_map.at(index)->first;
          const auto current_var =
              short_circuit_current_map.at((initial_index > a.first) ? std::make_pair(a.first, initial_index)
                                                                     : std::make_pair(initial_index, a.first));
          row[current_var] += (initial_index > a.first ? -1.0 : 1.0);
          continue;
        }

        row[index] += 1.0 / res;
        if (a.first != m_graph.begin()->first) {
          auto corrensponding_index = index_map.at(a.first);
          row[corrensponding_index] -= 1.0 / res;
        }

        row[size + num_short_circuits] -= emf / res;
      }
    }

    for (unsigned i = size; const auto &v : m_short_circuits) {
      auto row = extended_matrix[i++];

      if (v.first != m_graph.begin()->first) {
        auto first = index_map.at(v.first);
        row[first] = 1.0;
      }

      if (v.second != m_graph.begin()->first) {
        auto second = index_map.at(v.second);
        row[second] = -1.0;
      }

      row[size + num_short_circuits] = -v.emf;
    }

    return extended_matrix;
  };

  auto extended_matrix = make_extended_system();
  // Solve the linear system of equations to find unkown potentials and currents.
  auto unknowns = linmath::nonsingular_solver(linmath::matrix_d{std::move(extended_matrix)});

  auto result_potentials = solution_potentials{};
  // Fill base node potential with zero.
  result_potentials[m_graph.begin()->first] = 0.0;
  for (unsigned i = 0; i < size; ++i) {
    result_potentials[iterator_map[i]->first] = unknowns[i][0];
  }

  // Fill unkown currents that were found as a part of linear system of equations.
  auto result_currents = solution_currents{};
  for (const auto &v : short_circuit_current_map) {
    result_currents[v.first.first][v.first.second] = unknowns[v.second][0];
    result_currents[v.first.second][v.first.first] = -unknowns[v.second][0];
  }

  // Compute other currents from potentials, when there are no short-circuits.
  for (const auto &v : m_graph) {
    for (const auto &c : v.second) {
      if (throttle::is_roughly_equal(c.second.first, 0.0)) continue;
      result_currents[v.first][c.first] =
          (result_potentials[v.first] - result_potentials[c.first] + c.second.second) / c.second.first;
    }
  }

  return {result_potentials, result_currents};
}

} // namespace detail

void resistor_network::insert(unsigned first, unsigned second, double resistance, double emf) {
  resistance_emf_pair fwd_pair = {resistance, emf}, bck_pair = {resistance, -emf};
  m_graph.insert_edge({first, second}, fwd_pair, bck_pair);
}

std::vector<detail::connected_resistor_network> resistor_network::connected_components() const {
  auto components = m_graph.connected_components();
  return {components.begin(), components.end()};
}

solution resistor_network::solve() const {
  auto     components = connected_components();
  solution result;

  for (const auto &comp : components) {
    auto individual_sol = comp.solve();
    result.first.merge(individual_sol.first);
    result.second.merge(individual_sol.second);
  }

  return result;
}

} // namespace throttle::circuits