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
#include <boost/spirit/home/x3/support/utility/annotate_on_success.hpp>
#include <boost/spirit/home/x3/support/utility/error_reporting.hpp>

#include "linear_solver.hpp"
#include "resistor_network.hpp"

#include <range/v3/all.hpp>

namespace po = boost::program_options;
namespace x3 = boost::spirit::x3;
namespace ascii = x3::ascii;

namespace circuit_parser {

constexpr auto uint_ = x3::int_parser<unsigned>{};
using x3::double_;

namespace graph {
struct edge : x3::position_tagged {
  unsigned                first, second;
  double                  res;
  boost::optional<double> emf;
};
} // namespace graph

} // namespace circuit_parser

BOOST_FUSION_ADAPT_STRUCT(circuit_parser::graph::edge, (unsigned, first), (unsigned, second), (double, res),
                          (boost::optional<double>, emf))

namespace circuit_parser {

struct edge_class;
x3::rule<edge_class, circuit_parser::graph::edge> const edge = "edge";
const auto edge_def = uint_ > '-' > '-' > uint_ > ',' > double_ > ';' >> -(double_ >> 'V');

BOOST_SPIRIT_DEFINE(edge);

struct error_handler {
  x3::error_handler_result on_error(auto &first, auto const &last, auto const &x, auto const &context) {
    auto       &error_handler = x3::get<x3::error_handler_tag>(context).get();
    std::string message = "Error! Expecting: " + x.which() + " here:";
    error_handler(x.where(), message);
    return x3::error_handler_result::fail;
  }
};

struct edge_class : x3::annotate_on_success, error_handler {};

std::optional<std::pair<circuits::resistor_network, std::vector<circuit_parser::graph::edge>>> parse_circuit() {
  circuits::resistor_network               network;
  std::vector<circuit_parser::graph::edge> result;

  using ascii::space;

  using x3::phrase_parse;
  using x3::with;
  using iter_type = std::string::const_iterator;

  std::string input;
  std::noskipws(std::cin);
  std::copy(std::istream_iterator<char>{std::cin}, std::istream_iterator<char>{}, std::back_inserter(input));

  using error_handler_type = x3::error_handler<iter_type>;
  error_handler_type error_handler{input.begin(), input.end(), std::cerr};

  const auto edges = +edge;
  const auto parser = with<x3::error_handler_tag>(std::ref(error_handler))[edges];

  if (!phrase_parse(input.begin(), input.end(), parser, space, result)) return std::nullopt;

  for (const auto &v : result) {
    network.insert(v.first, v.second, v.res, v.emf.value_or(0.0));
  }

  return {std::make_pair(network, result)};
}

} // namespace circuit_parser

int main(int argc, char *argv[]) {
  auto parsed = circuit_parser::parse_circuit();

  if (!parsed) {
    std::cerr << "Aborting...\n";
    return EXIT_FAILURE;
  }

  auto [network, result] = parsed.value();
  auto [potentials, currents] = network.solve();

#if 0
  for (const auto &v : potentials) {
    std::cout << v.first << ": " << v.second << " V\n";
  }
#endif 

  for (const auto &v : result) {
    std::cout << v.first << " -- " << v.second << ": " << currents[v.first][v.second] << " A\n";
  }
}