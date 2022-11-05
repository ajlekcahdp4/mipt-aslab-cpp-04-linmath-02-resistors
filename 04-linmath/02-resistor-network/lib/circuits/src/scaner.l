/*
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <tsimmerman.ss@phystech.edu>, <alex.rom23@mail.ru> wrote this file.  As long as you
 * retain this notice you can do whatever you want with this stuff. If we meet
 * some day, and you think this stuff is worth it, you can buy me a beer in
 * return.
 * ----------------------------------------------------------------------------
 */

 %{
#include <iostream>

#include "bison_network_parser.hpp"
#include "scanner.hpp"

#define yyterminate() circuits::parser::make_EOF();

using namespace circuits;
%}

%option noyywrap nounput noinput nodefault 
%option c++
%option prefix="network_"
%option yyclass="scanner"


number      (0|[1-9][0-9]*)
float       [+-]?([0-9]*[.])?[0-9]+


%%

";"             { std::cout << "SEMICOL" << "\n"; return parser::make_SEMICOL(); }
","             { std::cout << "COLON" << "\n"; return parser::make_COMMA(); }
"--"            { std::cout << "LINE" << "\n"; return parser::make_LINE(); }
"V"             { std::cout << "VOLAGE" << "\n"; return parser::make_VOLTAGE(); }

{number}        { std::cout << yytext << "\n"; return parser::make_INTEGER(std::atoi(yytext)); }
{float}         { std::cout << yytext << "\n"; return parser::make_FLOAT(std::atof(yytext)); }

[ \t\n\r]       { }

<<EOF>>         { return yyterminate(); }

%%
