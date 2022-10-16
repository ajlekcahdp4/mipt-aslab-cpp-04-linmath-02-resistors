#include <boost/lexical_cast/bad_lexical_cast.hpp>
#include <cctype>
#include <chrono>
#include <cmath>
#include <iostream>
#include <iterator>
#include <set>
#include <string>

#include "contiguous_matrix.hpp"
#include "matrix.hpp"
#include "vector.hpp"

#include <concepts>
#include <optional>
#include <scn/scn.h>
#include <string>

#include <boost/fusion/adapted/std_tuple.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/option.hpp>
#include <boost/spirit/home/qi.hpp>
#include <boost/spirit/home/x3.hpp>

#include "linear_solver.hpp"
#include "resistor_network.hpp"

#include <range/v3/all.hpp>

namespace po = boost::program_options;
namespace x3 = boost::spirit::x3;
namespace ascii = x3::ascii;

std::optional<circuits::resistor_network> parse_circuit() {
  circuits::resistor_network network;

  constexpr auto uint_ = x3::int_parser<unsigned>{};

  using ascii::space;
  using x3::double_;
  using x3::phrase_parse;

  auto start = std::istream_iterator<char>{std::cin};
  auto end = std::istream_iterator<char>{};
  std::string input;
  std::copy(start, end, std::back_inserter(input));

  auto const edge_def = (uint_ >> '-' >> '-' >> uint_ >> ',' >> double_ >> ';' >> -(double_ >> 'V'));
  std::vector<std::tuple<unsigned, unsigned, double, boost::optional<double>>> result;

  bool r = phrase_parse(input.begin(), input.end(),
                        // Begin grammar
                        *edge_def,
                        // End grammar
                        space, result);

  if (!r) {
    std::cout << "Invalid input\n";
  }

  for (const auto &v : result) {
    network.insert(std::get<0>(v), std::get<1>(v), std::get<2>(v), std::get<3>(v).value_or(0.0));
  }

  return network;
}

int main(int argc, char *argv[]) {
  unsigned first, second;
  double res = 0, emf = 0;

  auto network = parse_circuit();
}