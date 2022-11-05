/*
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <tsimmerman.ss@phystech.edu>, <alex.rom23@mail.ru> wrote this file.  As long as you
 * retain this notice you can do whatever you want with this stuff. If we meet
 * some day, and you think this stuff is worth it, you can buy me a beer in
 * return.
 * ----------------------------------------------------------------------------
 */

%skeleton "lalr1.cc"
%require "3.8"

%defines

%define api.token.raw
%define api.parser.class { parser }
%define api.token.constructor
%define api.value.type variant
%define parse.assert
%define api.namespace { circuits }

%code requires {
#include <iostream>
#include <string>
#include <vector>

namespace circuits{
  class scanner;
  class driver;
}

using namespace circuits;

}

%code top
{

#include <iostream>
#include <string>

#include "frontend/scanner.hpp"
#include "bison_paracl_parser.hpp"
#include "frontend/driver.hpp"

static circuits::parser::symbol_type yylex(circuits::scanner &p_scanner, circuits::driver &p_driver) {
  return p_scanner.get_next_token();
}

}

%lex-param { circuits::scanner &scanner }
%lex-param { circuits::driver &driver }
%parse-param { circuits::scanner &scanner }
%parse-param { circuits::driver &driver }

%define parse.trace
%define parse.error verbose
%define api.token.prefix {TOKEN_}

/* Signle letter tokens */
%token LINE     "line"
%token SEMICOL  "semicol"
%token COMMA    "comma"

/* Terminals */
%token <int>    INTEGER "integer"
%token <float>  FLOAT   "float"

/* Rules that model the AST */
%type <int> vertex
%type <float> voltage

%start program

%%

program: statements { }


%%

// Bison expects us to provide implementation - otherwise linker complains
void circuits::parser::error(const std::string &message) {
  std::cout << "Error: " << message << "\n";
}