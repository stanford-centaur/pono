%{
/*********************                                                        */
/*! \file moxiparser.y
** \verbatim
** Top contributors (to current version):
**   Po-Chun Chien
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** (in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Bison grammar of the Model Exchange Interlingua (MoXI).
**
**        The actions hand what they recognize to a moxi::Script, which builds
**        the terms right away. Building them bottom-up, as the parser reduces,
**        keeps the depth of a term from ever reaching the call stack.
**
**/
%}

%skeleton "lalr1.cc"
%require "3.2"
%defines
%locations

%define api.namespace {pono::moxi}
%define api.parser.class {Parser}
%define api.location.type {pono::moxi::Location}
%define api.token.constructor
%define api.value.type variant
%define parse.error verbose

%param {pono::moxi::Script & drv}
%param {void * scanner}

%code requires {
  #include <string>
  #include <utility>
  #include <vector>

  #include "frontends/moxi_script.h"
}

%code {
  // The lexer, generated from moxiparser.l.
  pono::moxi::Parser::symbol_type moxilex(pono::moxi::Script & drv,
                                          void * scanner);
  #define yylex moxilex
}

%token END 0 "end of file"
%token LP "(" RP ")"

// reserved words
%token US "_" BANG "!" AS "as" LET "let" FORALL "forall" EXISTS "exists"

// commands
%token SET_LOGIC "set-logic" SET_INFO "set-info" SET_OPTION "set-option"
%token DECLARE_SORT "declare-sort" DEFINE_SORT "define-sort"
%token DECLARE_ENUM_SORT "declare-enum-sort"
%token DECLARE_CONST "declare-const" DEFINE_CONST "define-const"
%token DECLARE_FUN "declare-fun" DEFINE_FUN "define-fun"
%token DEFINE_SYSTEM "define-system" CHECK_SYSTEM "check-system" EXIT "exit"

// attributes of define-system and check-system
%token INPUT ":input" OUTPUT ":output" LOCAL ":local"
%token INIT ":init" TRANS ":trans" INV ":inv" SUBSYS ":subsys"
%token ASSUMPTION ":assumption" FAIRNESS ":fairness"
%token REACHABLE ":reachable" CURRENT ":current"
%token QUERY ":query" QUERIES ":queries"

%token <std::string> SYMBOL "symbol" PRIMED_SYMBOL "primed symbol"
%token <std::string> KEYWORD "keyword"
%token <std::string> NUMERAL "numeral" DECIMAL "decimal"
%token <std::string> BINARY "binary" HEXADECIMAL "hexadecimal"
%token <std::string> STRING "string"

%nterm <std::string> index
%nterm <std::vector<std::string>> symbol_list symbol_list1 index_list1
%nterm <pono::moxi::SortExpr> sort
%nterm <std::vector<pono::moxi::SortExpr>> sort_list sort_list1
%nterm <pono::moxi::SortedVar> sorted_var
%nterm <std::vector<pono::moxi::SortedVar>> sorted_var_list sorted_var_list1
%nterm <pono::moxi::Identifier> identifier qual_identifier
%nterm <smt::Term> term spec_constant
%nterm <smt::TermVec> term_list1
%nterm <pono::moxi::VarBinding> var_binding
%nterm <std::vector<pono::moxi::VarBinding>> var_binding_list1
%nterm <pono::moxi::QueryExpr> query
%nterm <std::vector<pono::moxi::QueryExpr>> query_list1

%%

script:
  %empty
| script command
;

command:
  "(" "set-logic" SYMBOL ")"
    { drv.set_logic($3, @3); }
| "(" "set-info" attribute ")"
| "(" "set-option" attribute ")"
| "(" "declare-sort" SYMBOL NUMERAL ")"
    { drv.declare_sort($3, $4, @3); }
| "(" "define-sort" SYMBOL "(" symbol_list ")" sort ")"
    { drv.define_sort($3, $5, $7, @3); }
| "(" "declare-enum-sort" SYMBOL "(" symbol_list1 ")" ")"
    { drv.declare_enum_sort($3, $5, @3); }
| "(" "declare-const" SYMBOL sort ")"
    { drv.declare_fun($3, {}, $4, @3); }
| "(" "declare-fun" SYMBOL "(" sort_list ")" sort ")"
    { drv.declare_fun($3, $5, $7, @3); }
| "(" "define-fun" SYMBOL "(" sorted_var_list ")" sort
    { drv.begin_define_fun($3, $5, $7, @3); }
  term ")"
    { drv.end_define_fun($9, @9); }
| "(" "define-const" SYMBOL sort
    { drv.begin_define_fun($3, {}, $4, @3); }
  term ")"
    { drv.end_define_fun($6, @6); }
| "(" "define-system" SYMBOL
    { drv.begin_system($3, @3); }
  variable_declarations
    { drv.end_variables(); }
  system_attributes ")"
    { drv.end_system(); }
| "(" "check-system" SYMBOL
    { drv.begin_check($3, @3); }
  variable_declarations
    { drv.end_variables(); }
  check_attributes ")"
    { drv.end_check(); }
| "(" "exit" ")"
    { YYACCEPT; }
| "(" SYMBOL
    { drv.error(@2, "unknown or unsupported command " + $2); }
  s_expr_list ")"
;

/* :input, :output and :local come before the other attributes, so that the
   formulas can refer to the variables they declare. */
variable_declarations:
  %empty
| variable_declarations variable_declaration
;

variable_declaration:
  ":input" "(" sorted_var_list ")"
    { drv.declare_variables(":input", $3, @1); }
| ":output" "(" sorted_var_list ")"
    { drv.declare_variables(":output", $3, @1); }
| ":local" "(" sorted_var_list ")"
    { drv.declare_variables(":local", $3, @1); }
;

system_attributes:
  %empty
| system_attributes system_attribute
;

system_attribute:
  ":init"
    { drv.begin_formula(":init"); }
  term
    { drv.set_system_formula(":init", $3, @1); }
| ":trans"
    { drv.begin_formula(":trans"); }
  term
    { drv.set_system_formula(":trans", $3, @1); }
| ":inv"
    { drv.begin_formula(":inv"); }
  term
    { drv.set_system_formula(":inv", $3, @1); }
| ":subsys" "(" SYMBOL "(" SYMBOL symbol_list ")" ")"
    { drv.add_subsystem($3, $5, $6, @3); }
;

check_attributes:
  %empty
| check_attributes check_attribute
;

check_attribute:
  ":assumption" "(" SYMBOL
    { drv.begin_formula(":assumption"); }
  term ")"
    { drv.add_check_formula(":assumption", $3, $5, @3); }
| ":fairness" "(" SYMBOL
    { drv.begin_formula(":fairness"); }
  term ")"
    { drv.add_check_formula(":fairness", $3, $5, @3); }
| ":reachable" "(" SYMBOL
    { drv.begin_formula(":reachable"); }
  term ")"
    { drv.add_check_formula(":reachable", $3, $5, @3); }
| ":current" "(" SYMBOL
    { drv.begin_formula(":current"); }
  term ")"
    { drv.add_check_formula(":current", $3, $5, @3); }
| ":query" query
    { drv.add_query($2); }
| ":queries" "(" query_list1 ")"
    {
      for (const auto & q : $3) {
        drv.add_query(q);
      }
    }
;

query:
  "(" SYMBOL "(" symbol_list ")" ")"
    { $$ = pono::moxi::QueryExpr{ $2, $4, @2 }; }
;

query_list1:
  query
    { $$.push_back(std::move($1)); }
| query_list1 query
    { $$ = std::move($1); $$.push_back(std::move($2)); }
;

sort:
  SYMBOL
    { $$ = pono::moxi::SortExpr{ $1, {}, {}, @1 }; }
| "(" "_" SYMBOL index_list1 ")"
    { $$ = pono::moxi::SortExpr{ $3, $4, {}, @$ }; }
| "(" SYMBOL sort_list1 ")"
    { $$ = pono::moxi::SortExpr{ $2, {}, $3, @$ }; }
;

sort_list:
  %empty
    {}
| sort_list sort
    { $$ = std::move($1); $$.push_back(std::move($2)); }
;

sort_list1:
  sort
    { $$.push_back(std::move($1)); }
| sort_list1 sort
    { $$ = std::move($1); $$.push_back(std::move($2)); }
;

sorted_var:
  "(" SYMBOL sort ")"
    { $$ = pono::moxi::SortedVar{ $2, $3, @2 }; }
;

sorted_var_list:
  %empty
    {}
| sorted_var_list sorted_var
    { $$ = std::move($1); $$.push_back(std::move($2)); }
;

sorted_var_list1:
  sorted_var
    { $$.push_back(std::move($1)); }
| sorted_var_list1 sorted_var
    { $$ = std::move($1); $$.push_back(std::move($2)); }
;

symbol_list:
  %empty
    {}
| symbol_list SYMBOL
    { $$ = std::move($1); $$.push_back(std::move($2)); }
;

symbol_list1:
  SYMBOL
    { $$.push_back(std::move($1)); }
| symbol_list1 SYMBOL
    { $$ = std::move($1); $$.push_back(std::move($2)); }
;

index:
  NUMERAL
    { $$ = $1; }
| SYMBOL
    { $$ = $1; }
;

index_list1:
  index
    { $$.push_back(std::move($1)); }
| index_list1 index
    { $$ = std::move($1); $$.push_back(std::move($2)); }
;

term:
  spec_constant
    { $$ = $1; }
| qual_identifier
    { $$ = drv.make_identifier_term($1); }
| "(" qual_identifier term_list1 ")"
    { $$ = drv.make_application($2, $3, @$); }
| "(" "let" "(" var_binding_list1 ")"
    { drv.push_let($4); }
  term ")"
    { drv.pop_scope(); $$ = $7; }
| "(" "forall" "(" sorted_var_list1 ")"
    { drv.push_quantifier($4); }
  term ")"
    { $$ = drv.pop_quantifier(true, $7, @$); }
| "(" "exists" "(" sorted_var_list1 ")"
    { drv.push_quantifier($4); }
  term ")"
    { $$ = drv.pop_quantifier(false, $7, @$); }
| "(" "!" term attribute_list1 ")"
    { $$ = $3; }
/* Not SMT-LIB, but files in the wild wrap terms in extra parentheses, which
   the MoXI tools written in Python accept too. */
| "(" term ")"
    { $$ = $2; }
;

term_list1:
  term
    { $$.push_back($1); }
| term_list1 term
    { $$ = std::move($1); $$.push_back($2); }
;

spec_constant:
  NUMERAL
    { $$ = drv.make_numeral($1, @1); }
| DECIMAL
    { $$ = drv.make_decimal($1, @1); }
| BINARY
    { $$ = drv.make_binary($1, @1); }
| HEXADECIMAL
    { $$ = drv.make_hexadecimal($1, @1); }
;

identifier:
  SYMBOL
    { $$ = pono::moxi::Identifier{ $1, {}, false, std::nullopt, @1 }; }
| PRIMED_SYMBOL
    { $$ = pono::moxi::Identifier{ $1, {}, true, std::nullopt, @1 }; }
| "(" "_" SYMBOL index_list1 ")"
    { $$ = pono::moxi::Identifier{ $3, $4, false, std::nullopt, @$ }; }
;

qual_identifier:
  identifier
    { $$ = std::move($1); }
| "(" "as" identifier sort ")"
    {
      $$ = std::move($3);
      $$.qualifier = std::move($4);
      $$.loc = @$;
    }
;

var_binding:
  "(" SYMBOL term ")"
    { $$ = pono::moxi::VarBinding{ $2, $3, @2 }; }
;

var_binding_list1:
  var_binding
    { $$.push_back(std::move($1)); }
| var_binding_list1 var_binding
    { $$ = std::move($1); $$.push_back(std::move($2)); }
;

/* Attributes of set-info, set-option and annotated terms are skipped. */
attribute:
  KEYWORD
| KEYWORD attribute_value
;

attribute_list1:
  attribute
| attribute_list1 attribute
;

attribute_value:
  atom
| "(" s_expr_list ")"
;

atom:
  NUMERAL
| DECIMAL
| BINARY
| HEXADECIMAL
| STRING
| SYMBOL
| PRIMED_SYMBOL
;

s_expr:
  atom
| KEYWORD
| "(" s_expr_list ")"
;

s_expr_list:
  %empty
| s_expr_list s_expr
;

%%

void pono::moxi::Parser::error(const location_type & loc,
                               const std::string & message)
{
  drv.error(loc, message);
}
