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
#include <boost/fusion/adapted/struct/adapt_struct.hpp>
#include <boost/fusion/include/adapt_struct.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>

#include <boost/program_options.hpp>
#include <boost/program_options/option.hpp>

#include <boost/spirit/home/x3.hpp>
#include <boost/spirit/home/x3/support/ast/position_tagged.hpp>

#include "linear_solver.hpp"
#include "resistor_network.hpp"

#include <range/v3/all.hpp>

namespace po = boost::program_options;
namespace x3 = boost::spirit::x3;
namespace ascii = x3::ascii;

namespace circuit_parser {

constexpr auto uint_ = x3::int_parser<unsigned>{};
constexpr auto double_ = x3::real_parser<double>{};

namespace graph {
struct edge : x3::position_tagged {
  unsigned first, second;
  double res;
  boost::optional<double> emf;
};
} // namespace graph

} // namespace circuit_parser

BOOST_FUSION_ADAPT_STRUCT(circuit_parser::graph::edge, (unsigned, first), (unsigned, second), (double, res),
                          (boost::optional<double>, emf))

namespace circuit_parser {

struct edge_class;
x3::rule<edge_class, circuit_parser::graph::edge> const edge = "edge";
const auto edge_def = uint_ > '-' > '-' > uint_ > ',' > double_ > ';' >> -(double_ > 'V');

BOOST_SPIRIT_DEFINE(edge);

std::optional<circuits::resistor_network> parse_circuit() {
  circuits::resistor_network network;
  std::vector<circuit_parser::graph::edge> result;
  
  using ascii::space;

  using x3::phrase_parse;
  using x3::with;
  using iter_type = std::string::const_iterator;

  std::string input;
  std::copy(std::istream_iterator<char>{std::cin}, std::istream_iterator<char>{}, std::back_inserter(input));

  bool r = phrase_parse(input.begin(), input.end(), +edge_def, space, result);

  if (!r) {
    std::cout << "Invalid input\n";
  }

  for (const auto &v : result) {
    network.insert(v.first, v.second, v.res, v.emf.value_or(0.0));
  }

  return network;
}

} // namespace circuit_parser

int main(int argc, char *argv[]) {
  unsigned first, second;
  double res = 0, emf = 0;

  auto network = circuit_parser::parse_circuit();
}