#include "linear_solver.hpp"
#include "matrix.hpp"
#include "vector.hpp"

#include <optional>
#include <utility>
#include <vector>

#include <map>
#include <unordered_map>

#pragma once

namespace circuits {

class resistor_network {
  using resistance_emf_pair = std::pair<double, double>;
  std::map<unsigned, std::map<unsigned, resistance_emf_pair>> m_map;

  struct short_circuit_edge {
    unsigned first;
    unsigned second;
    double   emf;
  };

  std::vector<short_circuit_edge> m_short_circuits;

public:
  void insert(unsigned first, unsigned second, double resistance, double emf);
  std::pair<std::unordered_map<unsigned, double>, std::unordered_map<unsigned, std::unordered_map<unsigned, double>>>
  solve() const;
};

} // namespace circuits