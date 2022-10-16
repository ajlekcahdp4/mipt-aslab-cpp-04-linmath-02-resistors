#include "resistor_network.hpp"
#include "matrix.hpp"

#include <stdexcept>
#include <utility>

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

}