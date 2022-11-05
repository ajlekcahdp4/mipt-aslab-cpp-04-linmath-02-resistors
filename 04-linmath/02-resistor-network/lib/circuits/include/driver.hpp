/*
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <tsimmerman.ss@phystech.edu>, <alex.rom23@mail.ru> wrote this file.  As long as you
 * retain this notice you can do whatever you want with this stuff. If we meet
 * some day, and you think this stuff is worth it, you can buy me a beer in
 * return.
 * ----------------------------------------------------------------------------
 */

#pragma once

#include "bison_network_parser.hpp"
#include "resistor_network.hpp"
#include "scanner.hpp"

#include <vector>

namespace circuits {

struct network_edge {
  unsigned                first, second;
  double                  res;
  boost::optional<double> emf;
};

class driver {
private:
  scanner m_scanner;
  parser  m_parser;

  friend class parser;
  friend class scanner;

public:
  std::pair<circuits::resistor_network, std::vector<network_edge>> m_parsed

  driver()
      : m_scanner{}, m_parser{m_scanner, *this} {}

  bool parse() { return m_parser.parse(); }
  void switch_input_stream(std::istream *is) { m_scanner.switch_streams(is, nullptr); }
};
} // namespace circuits