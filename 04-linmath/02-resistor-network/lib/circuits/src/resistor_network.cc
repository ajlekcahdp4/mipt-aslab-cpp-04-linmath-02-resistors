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

  std::unordered_map<unsigned, decltype(m_map)::const_iterator> iterator_map;
  std::unordered_map<unsigned, unsigned> index_map;

  unsigned j = 0;
  for (auto start = std::next(m_map.begin()), end = m_map.end(); start != end; ++start, ++j) {
    iterator_map[j] = start;
    index_map[start->first] = j;
  }

  auto size = iterator_map.size();

  throttle::linmath::contiguous_matrix<double> extended_matrix{size, size + 1};

  for (const auto &v : iterator_map) {
    const auto &[index, map_iter] = v;
    auto row = extended_matrix[index];

    for (const auto &a : map_iter->second) {
      auto [res, emf] = a.second;

      row[index] += 1.0 / res;
      if (a.first != m_map.begin()->first) {
        auto corrensponding_index = index_map[a.first];
        row[corrensponding_index] -= 1.0 / res;
      }
      
      row[size] -= emf / res;
    }
  }
  
  auto potentials = throttle::nonsingular_solver(std::move(extended_matrix));
  auto result = std::unordered_map<unsigned, double>{};

  result[m_map.begin()->first] = 0.0;
  for (unsigned i = 0; i < size; ++i) {
    result[iterator_map[i]->first] = potentials[i][0];
  }

  return result;
}

} // namespace circuits