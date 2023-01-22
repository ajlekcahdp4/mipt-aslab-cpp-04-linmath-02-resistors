#include "datastructures/vector.hpp"
#include "linmath/linear_solver.hpp"
#include "linmath/matrix.hpp"

#include <optional>
#include <utility>
#include <vector>

#include <unordered_map>
#include <unordered_set>

#pragma once

namespace circuits {

namespace linmath = throttle::linmath;

class connected_resistor_network {
  using resistance_emf_pair = std::pair<double, double>;
  std::unordered_map<unsigned, std::unordered_map<unsigned, resistance_emf_pair>> m_map;

  struct short_circuit_edge {
    unsigned first;
    unsigned second;
    double   emf;
  };

  std::vector<short_circuit_edge> m_short_circuits;

  void insert_impl(unsigned first, unsigned second, double resistance, double emf, bool to_throw);

public:
  void insert(unsigned first, unsigned second, double resistance, double emf);
  void try_insert(unsigned first, unsigned second, double resistance, double emf);

  using solution_potentials = std::unordered_map<unsigned, double>;
  using solution_currents = std::unordered_map<unsigned, std::unordered_map<unsigned, double>>;
  using solution = std::pair<solution_potentials, solution_currents>;

  solution solve() const;
};

class resistor_network {
public:
  using resistance_emf_pair = std::pair<double, double>;
  using map_type = std::unordered_map<unsigned, std::unordered_map<unsigned, resistance_emf_pair>>;

private:
  std::unordered_map<unsigned, std::unordered_map<unsigned, resistance_emf_pair>> m_map;

public:
  std::vector<connected_resistor_network> connected_components() const;
  const map_type                         &graph() const { return m_map; }

  void insert(unsigned first, unsigned second, double resistance = 0, double emf = 0);

  using solution_potentials = std::unordered_map<unsigned, double>;
  using solution_currents = std::unordered_map<unsigned, std::unordered_map<unsigned, double>>;
  using solution = std::pair<solution_potentials, solution_currents>;

  solution solve() const;
};

} // namespace circuits