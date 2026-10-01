%{
/*********************                                                        */
/*! \file moxiparser.l
** \verbatim
** Top contributors (to current version):
**   Po-Chun Chien
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** (in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Flex lexer of the Model Exchange Interlingua (MoXI), which follows
**        the lexical rules of SMT-LIB and adds the prime marking a
**        next-state variable, as in x' or |x|'.
**
**/

#include <cstdio>
#include <string>

#include "frontends/moxi_script.h"
#include "moxiparser.h"
#include "utils/exceptions.h"

// The parser passes its parameters on to the lexer; see moxiparser.y.
#define YY_DECL                                                  \
  pono::moxi::Parser::symbol_type moxilex(pono::moxi::Script & drv, \
                                          yyscan_t yyscanner)

// Every token extends the location by what it matched.
#define YY_USER_ACTION loc.columns(yyleng);

using pono::moxi::Parser;

namespace {

/** Moves the end of the location past the line breaks inside a token. */
void count_lines(pono::moxi::Location & loc, const char * text, int length)
{
  int lines = 0;
  int last_break = -1;
  for (int i = 0; i < length; ++i) {
    if (text[i] == '\n') {
      ++lines;
      last_break = i;
    }
  }
  if (lines) {
    loc.lines(lines);
    loc.columns(length - last_break - 1);
  }
}

/** @return the part of a token between its first and last few characters */
std::string strip(const char * text, int length, int front, int back)
{
  return std::string(text + front, static_cast<std::size_t>(length - front - back));
}

}  // namespace
%}

%option reentrant
%option noyywrap nounput noinput batch never-interactive
%option nodefault warn
%option prefix="moxi"

simple_start   [a-zA-Z~!@$%^&*_+=<>.?/\-]
simple_char    [a-zA-Z0-9~!@$%^&*_+=<>.?/\-]
simple_symbol  {simple_start}{simple_char}*
quoted_symbol  "|"[^|\\]*"|"
numeral        0|[1-9][0-9]*

%%

%{
  pono::moxi::Location & loc = drv.location();
  loc.step();
%}

[ \t\r\f\v]+          { loc.step(); }
\n+                   { loc.lines(yyleng); loc.step(); }
";"[^\n]*             { loc.step(); }

"("                   { return Parser::make_LP(loc); }
")"                   { return Parser::make_RP(loc); }

{numeral}             { return Parser::make_NUMERAL(yytext, loc); }
"0"[0-9]+             {
                        drv.error(loc, std::string("numeral with a leading zero: ")
                                           + yytext);
                      }
{numeral}"."[0-9]+    { return Parser::make_DECIMAL(yytext, loc); }
"#b"[01]+             { return Parser::make_BINARY(strip(yytext, yyleng, 2, 0), loc); }
"#x"[0-9a-fA-F]+      { return Parser::make_HEXADECIMAL(strip(yytext, yyleng, 2, 0), loc); }
\"([^"]|\"\")*\"      {
                        count_lines(loc, yytext, yyleng);
                        return Parser::make_STRING(strip(yytext, yyleng, 1, 1), loc);
                      }

":input"              { return Parser::make_INPUT(loc); }
":output"             { return Parser::make_OUTPUT(loc); }
":local"              { return Parser::make_LOCAL(loc); }
":init"               { return Parser::make_INIT(loc); }
":trans"              { return Parser::make_TRANS(loc); }
":inv"                { return Parser::make_INV(loc); }
":subsys"             { return Parser::make_SUBSYS(loc); }
":assumption"         { return Parser::make_ASSUMPTION(loc); }
":fairness"           { return Parser::make_FAIRNESS(loc); }
":reachable"          { return Parser::make_REACHABLE(loc); }
":current"            { return Parser::make_CURRENT(loc); }
":query"              { return Parser::make_QUERY(loc); }
":queries"            { return Parser::make_QUERIES(loc); }
":"{simple_symbol}    { return Parser::make_KEYWORD(yytext, loc); }

"_"                   { return Parser::make_US(loc); }
"!"                   { return Parser::make_BANG(loc); }
"as"                  { return Parser::make_AS(loc); }
"let"                 { return Parser::make_LET(loc); }
"forall"              { return Parser::make_FORALL(loc); }
"exists"              { return Parser::make_EXISTS(loc); }

"set-logic"           { return Parser::make_SET_LOGIC(loc); }
"set-info"            { return Parser::make_SET_INFO(loc); }
"set-option"          { return Parser::make_SET_OPTION(loc); }
"declare-sort"        { return Parser::make_DECLARE_SORT(loc); }
"define-sort"         { return Parser::make_DEFINE_SORT(loc); }
"declare-enum-sort"   { return Parser::make_DECLARE_ENUM_SORT(loc); }
"declare-const"       { return Parser::make_DECLARE_CONST(loc); }
"define-const"        { return Parser::make_DEFINE_CONST(loc); }
"declare-fun"         { return Parser::make_DECLARE_FUN(loc); }
"define-fun"          { return Parser::make_DEFINE_FUN(loc); }
"define-system"       { return Parser::make_DEFINE_SYSTEM(loc); }
"check-system"        { return Parser::make_CHECK_SYSTEM(loc); }
"exit"                { return Parser::make_EXIT(loc); }

{simple_symbol}"'"    { return Parser::make_PRIMED_SYMBOL(strip(yytext, yyleng, 0, 1), loc); }
{simple_symbol}       { return Parser::make_SYMBOL(yytext, loc); }
{quoted_symbol}"'"    {
                        count_lines(loc, yytext, yyleng);
                        return Parser::make_PRIMED_SYMBOL(strip(yytext, yyleng, 1, 2), loc);
                      }
{quoted_symbol}       {
                        count_lines(loc, yytext, yyleng);
                        return Parser::make_SYMBOL(strip(yytext, yyleng, 1, 1), loc);
                      }
"|"[^|\\]*            {
                        drv.error(loc, "quoted symbol without a closing |"
                                       " (a quoted symbol cannot contain \\)");
                      }

.                     {
                        drv.error(loc, std::string("unexpected character '")
                                           + yytext + "'");
                      }
<<EOF>>               { return Parser::make_END(loc); }

%%

void pono::moxi::Script::parse()
{
  // Releases the file and the lexer however parsing ends.
  struct Resources
  {
    std::FILE * file = nullptr;
    yyscan_t scanner = nullptr;
    ~Resources()
    {
      if (scanner) {
        moxilex_destroy(scanner);
      }
      if (file) {
        std::fclose(file);
      }
    }
  } resources;

  resources.file = std::fopen(filename_.c_str(), "r");
  if (!resources.file) {
    throw PonoException("Could not open MoXI file " + filename_);
  }
  if (moxilex_init(&resources.scanner)) {
    resources.scanner = nullptr;
    throw PonoException("Could not start the MoXI lexer");
  }
  moxiset_in(resources.file, resources.scanner);

  Parser parser(*this, resources.scanner);
  if (parser.parse()) {
    // The parser reports its errors by throwing, so this means it gave up.
    throw PonoException("Could not parse MoXI file " + filename_);
  }
}
