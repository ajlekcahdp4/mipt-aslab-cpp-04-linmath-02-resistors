#include "resistor_network.hpp"
#include "contiguous_matrix.hpp"
#include "disjoint_map_forest.hpp"
#include "equal.hpp"
#include "matrix.hpp"

#include <range/v3/all.hpp>

#include <stdexcept>
#include <utility>

#include <iterator>
#include <utility>

#include <unordered_map>
#include <unordered_set>

#include <boost/functional/hash.hpp>

namespace circuits {

void connected_resistor_network::insert(unsigned first, unsigned second, double resistance, double emf) {
  if (first == second) throw std::invalid_argument("Circuit graph can't have loops");

  auto found = m_map.find(first);
  if (found != m_map.end() && found->second.find(second) != found->second.end()) {
    throw std::invalid_argument("Edge is already present in the graph");
  }

  if (first > second) {
    std::swap(first, second);
    emf = -emf;
  }

  m_map[first].insert({second, std::make_pair(resistance, emf)});
  m_map[second].insert({first, std::make_pair(resistance, -emf)});

  if (throttle::is_roughly_equal(resistance, 0.0)) m_short_circuits.push_back({first, second, emf});
}

void connected_resistor_network::try_insert(unsigned first, unsigned second, double resistance, double emf) {
  if (first == second) throw std::invalid_argument("Circuit graph can't have loops");

  auto found = m_map.find(first);
  if (found != m_map.end() && found->second.find(second) != found->second.end()) {
    return;
  }

  if (first > second) {
    std::swap(first, second);
    emf = -emf;
  }

  m_map[first].insert({second, std::make_pair(resistance, emf)});
  m_map[second].insert({first, std::make_pair(resistance, -emf)});

  if (throttle::is_roughly_equal(resistance, 0.0)) m_short_circuits.push_back({first, second, emf});
}

connected_resistor_network::solution connected_resistor_network::solve() const {
  if (m_map.empty()) throw std::invalid_argument{"Network can't be empty"};

  // Maps indexes 0, 1, .... to corrensponding iterators in the unordered map that represents the input.
  std::unordered_map<unsigned, decltype(m_map)::const_iterator> iterator_map;
  // Maps indexes from input to 0, 1, ....
  std::unordered_map<unsigned, unsigned> index_map;
  // Maps pairs of input indexes to current variable
  std::unordered_map<std::pair<unsigned, unsigned>, unsigned, boost::hash<std::pair<unsigned, unsigned>>>
      short_circuit_current_map;

  unsigned j = 0;
  for (auto start = std::next(m_map.begin()), end = m_map.end(); start != end; ++start, ++j) {
    iterator_map[j] = start;
    index_map[start->first] = j;
  }

  const auto size = iterator_map.size();
  const auto num_short_circuits = m_short_circuits.size();

  for (j = size; const auto &v : m_short_circuits) {
    short_circuit_current_map.insert({{v.first, v.second}, j});
    ++j;
  }

  throttle::linmath::contiguous_matrix<double> extended_matrix{size + num_short_circuits,
                                                               size + num_short_circuits + 1};

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
      if (a.first != m_map.begin()->first) {
        auto corrensponding_index = index_map.at(a.first);
        row[corrensponding_index] -= 1.0 / res;
      }

      row[size + num_short_circuits] -= emf / res;
    }
  }

  for (unsigned i = size; const auto &v : m_short_circuits) {
    auto row = extended_matrix[i++];

    if (v.first != m_map.begin()->first) {
      auto first = index_map.at(v.first);
      row[first] = 1.0;
    }

    if (v.second != m_map.begin()->first) {
      auto second = index_map.at(v.second);
      row[second] = -1.0;
    }

    row[size + num_short_circuits] = -v.emf;
  }

#if 0
  for (unsigned i = 0; i < extended_matrix.rows(); ++i) {
    for (const auto &v : extended_matrix[i]) {
      std::cout << v << "\t";
    }
    std::cout << "\n";
  }
#endif

  // Solve the linear system of equations to find unkown potentials and currents.
  auto unknowns = throttle::nonsingular_solver(std::move(extended_matrix));

#if 0
  for (unsigned i = 0; i < unknowns.rows(); ++i) {
    for (const auto &v : unknowns[i]) {
      std::cout << v << "\t";
    }
    std::cout << "\n";
  }
#endif

  auto result_potentials = solution_potentials{};
  // Fill base node potential with zero.
  result_potentials[m_map.begin()->first] = 0.0;
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
  for (const auto &v : m_map) {
    for (const auto &c : v.second) {
      if (throttle::is_roughly_equal(c.second.first, 0.0)) continue;
      result_currents[v.first][c.first] =
          (result_potentials[v.first] - result_potentials[c.first] + c.second.second) / c.second.first;
    }
  }

  return {result_potentials, result_currents};
}

void resistor_network::insert(unsigned first, unsigned second, double resistance, double emf) {
  if (first == second) throw std::invalid_argument("Circuit graph can't have loops");

  if (first > second) {
    std::swap(first, second);
    emf = -emf;
  }

  auto found = m_map.find(first);
  if (found != m_map.end() && found->second.find(second) != found->second.end()) {
    throw std::invalid_argument("Edge is already present in the graph");
  }

  m_map[first].insert({second, std::make_pair(resistance, emf)});
  m_map[second].insert({first, std::make_pair(resistance, -emf)});
}

std::vector<connected_resistor_network> resistor_network::connected_components() const {
  throttle::disjoint_map_forest<unsigned, unsigned> dsu;

  for (const auto &v : m_map) {
    dsu.make_set(v.first, v.first);
  }

  for (const auto &v : m_map) {
    for (const auto &p : v.second) {
      // Here we do double work because each edge is visited twice, but this part isn't perfomance critical.
      dsu.union_set(v.first, p.first);
    }
  }

  std::unordered_map<unsigned, connected_resistor_network> connected_representatives;
  // FindSet does not change the component representaive. Here we iterate over all the nodes and find their
  // representaive to find all connected components.
  for (const auto &v : m_map) {
    auto found = dsu.find_set(v.first);
    if (connected_representatives.contains(*found)) continue;
    connected_representatives.insert({*found, connected_resistor_network{}});
  }

  for (const auto &v : m_map) {
    for (const auto &p : v.second) {
      auto found = dsu.find_set(v.first);
      connected_representatives[*found].try_insert(v.first, p.first, p.second.first, p.second.second);
    }
  }

  std::vector<connected_resistor_network> result;
  for (const auto &v : connected_representatives) {
    result.push_back(std::move(v.second));
  }

  return result;
}

resistor_network::solution resistor_network::solve() const {
  auto     components = connected_components();
  solution result;

  for (const auto &comp: components) {
    auto individual_sol = comp.solve();
    result.first.merge(individual_sol.first);
    result.second.merge(individual_sol.second);
  }

  return result;
}

} // namespace circuits