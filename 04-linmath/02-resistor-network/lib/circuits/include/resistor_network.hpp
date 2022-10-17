#include "linear_solver.hpp"
#include "matrix.hpp"
#include "vector.hpp"

#include <utility>
#include <optional>
#include <vector>

#include <map>
#include <unordered_map>

#pragma once

namespace circuits {

class resistor_network {
  using resistance_emf_pair = std::pair<double, double>;
  std::map<unsigned, std::map<unsigned, resistance_emf_pair>> m_map;

public:
  void insert(unsigned first, unsigned second, double resistance, double emf);
  std::unordered_map<unsigned, double> solve() const;
};

} // namespace circuits