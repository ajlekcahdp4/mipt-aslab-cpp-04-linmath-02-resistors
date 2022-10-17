#include "resistor_network.hpp"
#include "contiguous_matrix.hpp"
#include "equal.hpp"
#include "matrix.hpp"

#include <range/v3/all.hpp>

#include <stdexcept>
#include <utility>

#include <iterator>
#include <utility>

namespace circuits {

void resistor_network::insert(unsigned int first, unsigned int second, double resistance, double emf) {
  if (first == second) throw std::invalid_argument("Circuit graph can't have loops");

  auto found = m_map.find(first);
  if (found != m_map.end() && found->second.find(second) != found->second.end()) {
    throw std::invalid_argument("Edge is already present in the graph");
  }

  m_map[first].insert({second, std::make_pair(resistance, emf)});
  m_map[second].insert({first, std::make_pair(resistance, -emf)});
}

std::unordered_map<unsigned, double> resistor_network::solve() const {
  if (m_map.empty()) throw std::invalid_argument{"Network can't be empty"};

  unsigned min_index = m_map.begin()->first;
  unsigned max_index = std::prev(m_map.end())->first;

  auto                                         size = max_index - min_index + 1;
  throttle::linmath::contiguous_matrix<double> extended_matrix{size, size + 1};

  for (const auto &vertex : m_map) {
    const auto  index = vertex.first - min_index;
    const auto &map = vertex.second;

    for (const auto &another : map) {
      const auto another_index = another.first - min_index;
      const auto [res, emf] = another.second;
      auto row = extended_matrix[index];

      if (throttle::is_roughly_equal(res, 0.0)) {
        ranges::fill(row, 0);

        row[index] += (index == 0 ? 0.0 : 1.0);
        row[another_index] -= 1.0;
        row[size] -= emf;

        break;
      }

      row[index] += (index == 0 ? 0.0 : 1.0 / res);
      row[another_index] -= 1.0 / res;
      row[size] -= emf / res;
    }
  }

  for (unsigned i = 0; i < size; ++i) {
    for (const auto &v : extended_matrix[i]) {
      std::cout << v << "\t";
    }
    std::cout << "\n";
  }

  auto potentials = throttle::nonsingular_solver(std::move(extended_matrix));
  auto result = std::unordered_map<unsigned, double>{};

  for (unsigned i = 0; i < size; ++i) {
    result[i + min_index] = potentials[i][0];
  }

  return result;
}

} // namespace circuits