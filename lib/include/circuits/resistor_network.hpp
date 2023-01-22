/*
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <tsimmerman.ss@phystech.edu>, <alex.rom23@mail.ru> wrote this file.  As long as you
 * retain this notice you can do whatever you want with this stuff. If we meet
 * some day, and you think this stuff is worth it, you can buy us a beer in
 * return.
 * ----------------------------------------------------------------------------
 */

#include "datastructures/ud_asymmetric_graph.hpp"
#include "datastructures/vector.hpp"
#include "linmath/linear_solver.hpp"
#include "linmath/matrix.hpp"

#include <optional>
#include <utility>
#include <vector>

#include <unordered_map>
#include <unordered_set>

#pragma once

namespace throttle::circuits {

using resistance_emf_pair = std::pair<double, double>;

struct short_circuit_edge {
  unsigned first;
  unsigned second;
  double   emf;
};

using circuit_graph_type = containers::ud_asymmetric_graph<unsigned, resistance_emf_pair>;
using solution_potentials = std::unordered_map<unsigned, double>;
using solution_currents = std::unordered_map<unsigned, std::unordered_map<unsigned, double>>;
using solution = std::pair<solution_potentials, solution_currents>;

namespace detail {
class connected_resistor_network {
  circuit_graph_type              m_graph;
  std::vector<short_circuit_edge> m_short_circuits;

public:
  connected_resistor_network(circuit_graph_type graph);

  solution solve() const;
};
} // namespace detail

class resistor_network {
  circuit_graph_type m_graph;

public:
  std::vector<detail::connected_resistor_network> connected_components() const;

  circuit_graph_type graph() const { return m_graph; }
  solution           solve() const;

  void insert(unsigned first, unsigned second, double resistance = 0, double emf = 0);
};

} // namespace throttle::circuits