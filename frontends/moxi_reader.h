/*********************                                                        */
/*! \file moxi_reader.h
** \verbatim
** Top contributors (to current version):
**   Po-Chun Chien
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** (in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Reads a MoXI file into the systems and checks it defines, as
**        terms of a solver, for MoxiEncoder to encode.
**
**        The Flex/Bison parser hands what it recognizes to a moxi::Reader,
**        which resolves the names, checks the sorts and builds the terms in
**        the solver right away. The variables of a system stand in its terms
**        as placeholder symbols, which MoxiEncoder substitutes when it
**        instantiates the system.
**
**/

#pragma once

#include <cstddef>
#include <cstdint>
#include <deque>
#include <iosfwd>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "smt-switch/smt.h"

namespace pono {
namespace moxi {

/** A point in the input, counted from 1. */
struct Position
{
  std::uint32_t line = 1;
  std::uint32_t column = 1;
};

/** A span of the input. It serves as the location type of the parser, which
 *  expects the begin and end members. */
struct Location
{
  Position begin;
  Position end;

  /** Starts the next span where this one ends. */
  void step() { begin = end; }

  /** Extends the span by some columns on its last line. */
  void columns(int count);

  /** Extends the span by some lines, to the start of the last of them. */
  void lines(int count);
};

std::ostream & operator<<(std::ostream & os, const Location & loc);

/** A sort as written, before its names are resolved. */
struct SortExpr
{
  std::string name;
  std::vector<std::string> indices;  ///< e.g. the width in (_ BitVec 8)
  std::vector<SortExpr> args;        ///< e.g. the sorts in (Array Int Bool)
  Location loc;
};

// The parser holds its values in slots as large as the largest of their
// types, clears each slot it makes, and moves the values between slots.
// Hence the sorts nested in the values below are kept behind pointers.

/** An identifier as written: a symbol, possibly indexed, primed, or
 *  qualified with a sort as in (as const (Array Int Int)). */
struct Identifier
{
  std::string name;
  std::vector<std::string> indices;
  bool primed = false;
  std::shared_ptr<const SortExpr> qualifier;  ///< null if unqualified
  Location loc;
};

/** A symbol declared with a sort, e.g. a variable or a parameter. */
struct SortedVar
{
  std::string name;
  std::shared_ptr<const SortExpr> sort;  ///< never null
  Location loc;
};

/** A symbol that a let binds to a term. */
struct VarBinding
{
  std::string name;
  smt::Term term;
  Location loc;
};

/** A query of a check-system command, before its names are resolved. */
struct QueryExpr
{
  std::string name;
  std::vector<std::string> formulas;
  Location loc;
};

/** An enumeration sort. Its values are encoded as the bit-vectors 0, 1, ...
 *  of the least width that fits them all. */
struct EnumSort
{
  std::string name;
  std::vector<std::string> values;
  smt::Sort sort;
};

/** A sort resolved against the declarations in scope. The solver sort alone
 *  cannot tell an enumeration from the bit-vectors encoding it. */
struct SortInfo
{
  smt::Sort sort;
  const EnumSort * enumeration = nullptr;  ///< set if the sort is an enum

  bool operator==(const SortInfo & other) const
  {
    return sort == other.sort && enumeration == other.enumeration;
  }
  bool operator!=(const SortInfo & other) const { return !(*this == other); }
};

/** A variable of a system or check. Its current and next values are
 *  placeholder symbols in the terms of the command that declares it. */
struct Variable
{
  std::string name;
  SortInfo sort;
  smt::Term curr;
  smt::Term next;
  Location loc;
};

/** A symbol declared with declare-const, or declare-fun without arguments.
 *  Its value is fixed over time, and it has a placeholder symbol like a
 *  variable does. */
struct Constant
{
  std::string name;
  SortInfo sort;
  smt::Term placeholder;
  Location loc;
};

struct System;

/** An instance of a system inside another, from a :subsys attribute. */
struct Subsystem
{
  std::string name;
  std::shared_ptr<const System> system;
  /** For each input and then output of the instantiated system, the position
   *  of the argument among the variables of the enclosing system, in the
   *  order System::variable counts them. */
  std::vector<std::size_t> args;
  Location loc;
};

/** A system, from a define-system command. Absent formulas are null. */
struct System
{
  std::string name;
  std::vector<Variable> inputs;
  std::vector<Variable> outputs;
  std::vector<Variable> locals;
  smt::Term init;
  smt::Term trans;
  smt::Term inv;
  std::vector<Subsystem> subsystems;
  Location loc;

  std::size_t num_variables() const
  {
    return inputs.size() + outputs.size() + locals.size();
  }

  /** @return the variable at position i, counting inputs, then outputs, then
   *  locals */
  const Variable & variable(std::size_t i) const;
};

enum class FormulaKind
{
  ASSUMPTION,
  FAIRNESS,
  REACHABLE,
  CURRENT
};

/** A named formula of a check-system command, e.g. from :reachable. */
struct CheckFormula
{
  std::string name;
  FormulaKind kind;
  smt::Term term;
  Location loc;
};

/** A query of a check-system command, from :query or one of :queries. */
struct Query
{
  std::string name;
  std::vector<std::size_t> formulas;  ///< positions in Check::formulas
  Location loc;
};

/** A check-system command. Its variables rename, in order, the inputs,
 *  outputs and locals of the system it checks. */
struct Check
{
  std::shared_ptr<const System> system;
  std::vector<Variable> inputs;
  std::vector<Variable> outputs;
  std::vector<Variable> locals;
  std::vector<CheckFormula> formulas;
  std::vector<Query> queries;
  Location loc;
};

/** Reads a MoXI file: parses it and holds what it defines.
 *
 *  Besides the results, it offers the operations that the parser and the
 *  lexer build them with. Errors are reported by throwing a PonoException
 *  that names the place in the file.
 */
class Reader
{
 public:
  /** Parses a MoXI file.
   *  @param filename the file to read
   *  @param solver the solver to build the terms and sorts with
   */
  Reader(const std::string & filename, const smt::SmtSolver & solver);

  Reader(const Reader &) = delete;
  Reader & operator=(const Reader &) = delete;

  const std::string & filename() const { return filename_; }

  const smt::SmtSolver & solver() const { return solver_; }

  /** @return the check-system commands, in the order of the file */
  const std::vector<Check> & checks() const { return checks_; }

  /** @return the declared constants, in the order of the file */
  const std::deque<Constant> & constants() const { return constants_; }

  /** @return the names of the uninterpreted functions, which are symbols of
   *  the solver themselves */
  const std::unordered_set<std::string> & function_names() const
  {
    return function_names_;
  }

  /** Throws a PonoException locating the message in the file. */
  [[noreturn]] void error(const Location & loc,
                          const std::string & message) const;

  // The rest is for the lexer and the parser.

  /** @return the location the lexer keeps up to date */
  Location & location() { return location_; }

  // commands
  void set_logic(const std::string & logic, const Location & loc);
  void declare_sort(const std::string & name,
                    const std::string & arity,
                    const Location & loc);
  void define_sort(const std::string & name,
                   const std::vector<std::string> & params,
                   const SortExpr & body,
                   const Location & loc);
  void declare_enum_sort(const std::string & name,
                         const std::vector<std::string> & values,
                         const Location & loc);
  void declare_fun(const std::string & name,
                   const std::vector<SortExpr> & args,
                   const SortExpr & result,
                   const Location & loc);
  void begin_define_fun(const std::string & name,
                        const std::vector<SortedVar> & params,
                        const SortExpr & result,
                        const Location & loc);
  void end_define_fun(const smt::Term & body, const Location & loc);

  void begin_system(const std::string & name, const Location & loc);
  void end_system();
  void begin_check(const std::string & system, const Location & loc);
  void end_check();

  /** Declares the variables of an :input, :output or :local attribute, of
   *  the system or check being parsed. */
  void declare_variables(const std::string & attribute,
                         const std::vector<SortedVar> & vars,
                         const Location & loc);
  /** Marks the end of the variable declarations, which come first. */
  void end_variables();

  /** Prepares for the formula of an attribute, e.g. :init, which decides
   *  whether next-state variables may appear in it. */
  void begin_formula(const std::string & attribute);
  void set_system_formula(const std::string & attribute,
                          const smt::Term & formula,
                          const Location & loc);
  void add_subsystem(const std::string & name,
                     const std::string & system,
                     const std::vector<std::string> & args,
                     const Location & loc);
  void add_check_formula(const std::string & attribute,
                         const std::string & name,
                         const smt::Term & formula,
                         const Location & loc);
  void add_query(const QueryExpr & query);

  // terms
  smt::Term make_numeral(const std::string & digits, const Location & loc);
  smt::Term make_decimal(const std::string & text, const Location & loc);
  smt::Term make_binary(const std::string & bits, const Location & loc);
  smt::Term make_hexadecimal(const std::string & digits, const Location & loc);
  /** @return the term an identifier stands for on its own */
  smt::Term make_identifier_term(const Identifier & id);
  /** @return the application of the function an identifier names */
  smt::Term make_application(const Identifier & id,
                             const smt::TermVec & args,
                             const Location & loc);
  /** Binds the symbols of a let for its body. */
  void push_let(const std::vector<VarBinding> & bindings);
  /** Ends the scope of the innermost let, or of define-fun parameters. */
  void pop_scope();
  /** Binds the variables of a quantifier for its body. */
  void push_quantifier(const std::vector<SortedVar> & vars);
  /** Ends the scope of the innermost quantifier.
   *  @return the quantified formula */
  smt::Term pop_quantifier(bool forall,
                           const smt::Term & body,
                           const Location & loc);

 private:
  /** Reads the file with the lexer and parser, defined in moxiparser.l. */
  void parse();

  /** @return a symbol of the solver to stand for a value in the terms of the
   *  file, named apart from any other that any reader makes */
  smt::Term make_placeholder(const smt::Sort & sort);

  /** A function from define-fun, applied by substituting its body. */
  struct Macro
  {
    smt::TermVec params;
    std::vector<SortInfo> param_sorts;
    smt::Term body;
    SortInfo result;
  };

  /** An uninterpreted function from declare-fun. */
  struct Function
  {
    smt::Term symbol;
    std::vector<SortInfo> arg_sorts;
    SortInfo result;
  };

  /** A sort from define-sort, possibly with parameters. */
  struct SortDefinition
  {
    std::vector<std::string> params;
    SortExpr body;
  };

  /** A sort from declare-sort: the sort, or the constructor if the arity is
   *  positive. */
  struct DeclaredSort
  {
    smt::Sort sort;
    std::uint64_t arity;
  };

  /** What the command being parsed is, for the attributes it may take. */
  enum class Scope
  {
    TOP,
    SYSTEM,
    CHECK
  };

  using SortParams = std::unordered_map<std::string, SortInfo>;

  SortInfo resolve_sort(const SortExpr & expr,
                        const SortParams * params = nullptr);
  smt::Sort int_sort(const Location & loc);
  smt::Sort real_sort(const Location & loc);
  /** @return the sort that numerals have in the logic */
  smt::Sort numeral_sort(const Location & loc);
  std::uint64_t parse_index(const std::string & index,
                            const Location & loc) const;

  /** Refuses a name for a new sort that is taken. */
  void check_new_sort_name(const std::string & name,
                           const Location & loc) const;
  /** Refuses a name for a new global symbol that is taken. */
  void check_new_global_name(const std::string & name,
                             const Location & loc) const;

  /** @return the variable of the system or check in scope, or null */
  const Variable * lookup_variable(const std::string & name) const;
  /** @return the term a let or quantifier binds the name to, or null */
  smt::Term lookup_bound(const std::string & name) const;
  void bind(const std::string & name, const smt::Term & term);
  smt::Term make_next(const std::string & name, const Location & loc);

  smt::Term make_indexed_constant(const Identifier & id);
  smt::Term make_indexed_application(const Identifier & id,
                                     const smt::TermVec & args,
                                     const Location & loc);
  smt::Term make_builtin(const std::string & name,
                         smt::TermVec args,
                         const Location & loc);
  smt::Term apply_macro(const std::string & name,
                        const Macro & macro,
                        smt::TermVec args,
                        const Location & loc);
  /** Converts the integer arguments to reals if any argument is a real. */
  void unify_numeric(smt::TermVec & args, const Location & loc) const;
  /** Converts an integer term to a real one if the sort asks for that. */
  smt::Term coerce(const smt::Term & term,
                   const smt::Sort & sort,
                   const Location & loc) const;
  /** Builds a term, reporting the solver's objections at the location. */
  smt::Term make_term(const smt::Op & op,
                      const smt::TermVec & args,
                      const Location & loc) const;

  /** @return the variables declared so far by the command being parsed */
  std::vector<Variable> & variables(const std::string & attribute,
                                    const Location & loc);
  /** Checks that the variables a check-system declares rename those of the
   *  system, or declares them with the system's names if it declares none. */
  void check_signature(const std::string & attribute,
                       std::vector<Variable> & declared,
                       const std::vector<Variable> & expected,
                       const Location & loc);
  /** Makes the variables of the command being parsed visible to its terms. */
  void enter_variables(const std::vector<Variable> & vars);

  std::string filename_;
  smt::SmtSolver solver_;
  Location location_;

  // the logic
  bool logic_set_ = false;
  bool logic_has_int_ = true;
  bool logic_has_real_ = true;
  smt::Sort bool_sort_;
  smt::Sort int_sort_;
  smt::Sort real_sort_;

  // global declarations
  std::unordered_map<std::string, DeclaredSort> declared_sorts_;
  std::unordered_map<std::string, SortDefinition> defined_sorts_;
  std::deque<EnumSort> enum_sorts_;  ///< a deque keeps pointers valid
  std::unordered_map<std::string, const EnumSort *> enum_sort_names_;
  std::unordered_map<std::string, smt::Term> enum_values_;
  std::deque<Constant> constants_;
  std::unordered_map<std::string, const Constant *> constant_names_;
  std::unordered_map<std::string, Function> functions_;
  std::unordered_set<std::string> function_names_;
  std::unordered_map<std::string, Macro> macros_;
  std::unordered_map<std::string, std::shared_ptr<const System>> systems_;
  std::vector<Check> checks_;

  // the command being parsed
  Scope scope_ = Scope::TOP;
  std::shared_ptr<System> system_;  ///< while parsing a define-system
  Check check_;                     ///< while parsing a check-system
  std::unordered_set<std::string> declared_attributes_;
  std::unordered_set<std::string> command_names_;  ///< its names in use
  std::unordered_map<std::string, const Variable *> variables_in_scope_;
  /** positions of the system's variables, as System::variable counts them */
  std::unordered_map<std::string, std::size_t> variable_positions_;
  /** positions of the check's formulas in Check::formulas */
  std::unordered_map<std::string, std::size_t> formula_positions_;
  std::vector<QueryExpr> query_exprs_;  ///< the check's queries, unresolved
  bool next_allowed_ = false;           ///< whether the formula may use primes
  std::string formula_attribute_;
  std::string macro_name_;  ///< of the define-fun being parsed
  Macro macro_;

  // symbols bound by lets and quantifiers, innermost last
  std::unordered_map<std::string, std::vector<smt::Term>> bound_;
  std::vector<std::vector<std::string>> bound_scopes_;
  std::vector<smt::TermVec> quantifier_params_;
};

}  // namespace moxi
}  // namespace pono
