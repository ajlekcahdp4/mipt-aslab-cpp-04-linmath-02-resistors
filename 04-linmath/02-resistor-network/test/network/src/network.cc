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

#include <boost/fusion/adapted/std_pair.hpp>
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

namespace graph {
struct network_edge : x3::position_tagged {
  unsigned                first, second;
  double                  res;
  boost::optional<double> emf;
};

} // namespace graph

struct error_handler {
  x3::error_handler_result on_error(auto &, auto const &, auto const &x, auto const &context) {
    auto       &error_handler = x3::get<x3::error_handler_tag>(context).get();
    std::string message = "Parsing error! Expecting: '" + x.which() + "' here:";
    error_handler(x.where(), message);
    return x3::error_handler_result::fail;
  }
};

} // namespace circuit_parser

BOOST_FUSION_ADAPT_STRUCT(circuit_parser::graph::network_edge, first, second, res, emf)

namespace circuit_parser {

struct rule_d : error_handler {};
struct rule_u : error_handler {};

constexpr x3::rule<rule_d, double> double_named = {"double"}; 
constexpr x3::rule<rule_u, unsigned> unsigned_named = {"unsigned"};

constexpr auto double_named_def = x3::real_parser<double>{};
constexpr auto unsigned_named_def = x3::int_parser<unsigned>{};

BOOST_SPIRIT_DEFINE(double_named, unsigned_named);

struct edge_class;

constexpr x3::rule<edge_class, circuit_parser::graph::network_edge> const edge = "edge";

const auto edge_def = unsigned_named > '-' > '-' > unsigned_named > ',' > double_named > ';' >> -(double_named >> 'V');

BOOST_SPIRIT_DEFINE(edge);

struct edge_class : x3::annotate_on_success, error_handler{};

std::optional<std::pair<circuits::resistor_network, std::vector<circuit_parser::graph::network_edge>>> parse_circuit() {
  std::vector<circuit_parser::graph::network_edge> parse_result;

  using ascii::space;

  using x3::phrase_parse;
  using x3::with;
  using iter_type = std::string::const_iterator;

  std::string input;
  std::noskipws(std::cin);
  std::copy(std::istream_iterator<char>{std::cin}, std::istream_iterator<char>{}, std::back_inserter(input));

  using error_handler_type = x3::error_handler<iter_type>;
  error_handler_type error_handler{input.begin(), input.end(), std::cerr};

  const auto parser = with<x3::error_handler_tag>(std::ref(error_handler))[+edge];
  bool       res = phrase_parse(input.begin(), input.end(), parser, space, parse_result);

  if (!res) return std::nullopt;

  circuits::resistor_network network;

  for (const auto &v : parse_result) {
    network.insert(v.first, v.second, v.res, v.emf.value_or(0.0));
  }

  return std::make_pair(std::move(network), std::move(parse_result));
}

} // namespace circuit_parser

int main(int argc, char *argv[]) {
  bool non_verbose = false, pot_verbose = false;

  po::options_description desc("Available options");
  desc.add_options()("help,h", "Print this help message")("nonverbose,n", "Non-verbose output")("potentials,p",
                                                                                         "Print vertex potentials");
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  po::notify(vm);

  non_verbose = vm.count("nonverbose");
  if (vm.count("potentials")) {
    pot_verbose = true;
    non_verbose = false;
  }

  auto parsed = circuit_parser::parse_circuit();

  if (!parsed) {
    std::cerr << "Aborting...\n";
    return EXIT_FAILURE;
  }

  auto [network, result] = parsed.value();
  auto [potentials, currents] = network.solve();

  for (const auto &v : result) {
    if (!non_verbose) {
      std::cout << v.first << " -- " << v.second << ": " << currents[v.first][v.second] << " A\n";
    } else {
      std::cout << currents[v.first][v.second] << "\n";
    }
  }

  if (pot_verbose) {
    for (const auto &v : result) {
      std::cout << v.first << " -- " << potentials[v.first] << " V\n";
    }
  }
}