#include "linear_solver.hpp"
#include "matrix.hpp"
#include "vector.hpp"

#include <optional>
#include <utility>
#include <vector>

#include <unordered_set>
#include <unordered_map>

#pragma once

namespace circuits {

class connected_resistor_network {
  using resistance_emf_pair = std::pair<double, double>;
  std::unordered_map<unsigned, std::map<unsigned, resistance_emf_pair>> m_map;

  struct short_circuit_edge {
    unsigned first;
    unsigned second;
    double   emf;
  };

  std::vector<short_circuit_edge> m_short_circuits;

public:
  void insert(unsigned first, unsigned second, double resistance, double emf);
  void try_insert(unsigned first, unsigned second, double resistance, double emf);

  using solution_potentials = std::unordered_map<unsigned, double>;
  using solution_currents = std::unordered_map<unsigned, std::unordered_map<unsigned, double>>;
  using solution = std::pair<solution_potentials, solution_currents>;

  solution solve() const;
};

class resistor_network {
  using resistance_emf_pair = std::pair<double, double>;
  std::unordered_map<unsigned, std::unordered_map<unsigned, resistance_emf_pair>> m_map;

public:
  std::vector<connected_resistor_network> connected_components() const;
  void insert(unsigned first, unsigned second, double resistance, double emf);

  using solution_potentials = std::unordered_map<unsigned, double>;
  using solution_currents = std::unordered_map<unsigned, std::unordered_map<unsigned, double>>;
  using solution = std::pair<solution_potentials, solution_currents>;

  solution solve() const;
};

} // namespace circuits