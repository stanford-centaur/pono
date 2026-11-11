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

#include "options/options.h"
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

/** An enumeration sort that another sort encodes, its values being 0, 1, ...
 *  in the order of the declaration. */
struct Enumeration
{
  std::string name;
  smt::Sort sort;  ///< that encodes it
  smt::TermVec values;
  /** whether the encoding sort has no other values, as the bit-vectors of a
   *  width do if there are 2^width values */
  bool exact = false;
};

struct EncodedSort;

/** A sort of MoXI, as the solver's sort that encodes it. An enumeration that
 *  another sort encodes, e.g. Int, has to be told apart from that sort and
 *  from other enumerations, and so do the sorts made of it. */
struct MoxiSort
{
  smt::Sort sort;
  /** what tells the sort apart if it involves such an enumeration; null
   *  otherwise, as the solver's sort is then the same */
  const EncodedSort * encoded = nullptr;

  /** @return whether the sort involves an enumeration that another sort
   *  encodes */
  bool involves_enumeration() const { return encoded; }

  std::string to_string() const;

  bool operator==(const MoxiSort & other) const
  {
    return sort == other.sort && encoded == other.encoded;
  }

  bool operator!=(const MoxiSort & other) const { return !(*this == other); }
};

/** A sort that involves an enumeration that another sort encodes. The reader
 *  holds one of each, so that two are the same sort exactly if they are the
 *  same object. */
struct EncodedSort
{
  smt::Sort sort;    ///< that encodes it
  std::string name;  ///< as MoXI writes it
  /** the enumeration the sort is, or null */
  const Enumeration * enumeration = nullptr;
  /** the index and element sorts of an array, or the arguments of a declared
   *  sort */
  std::vector<MoxiSort> args;
};

/** A term of the solver, with its sort of MoXI. */
struct MoxiTerm
{
  smt::Term term;
  /** what tells its sort apart if it involves an enumeration that another
   *  sort encodes; null otherwise, as the term's sort is then its sort of
   *  MoXI. Null, rather than an optional sort, keeps the term as cheap to
   *  copy as the solver's own term is, as nearly all terms need nothing. */
  const EncodedSort * encoded = nullptr;

  /** @return its sort of MoXI */
  MoxiSort sort() const;
};

using MoxiTermVec = std::vector<MoxiTerm>;

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
  MoxiTerm term;
  Location loc;
};

/** A query of a check-system command, before its names are resolved. */
struct QueryExpr
{
  std::string name;
  std::vector<std::string> formulas;
  Location loc;
};

/** A variable of a system or check. Its current and next values are
 *  placeholder symbols in the terms of the command that declares it. */
struct Variable
{
  std::string name;
  MoxiSort sort;
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
  MoxiSort sort;
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
  std::vector<Location> queries_attributes;  ///< where its :queries are
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
   *  @param enum_encoding the sort to encode enumeration sorts as
   */
  Reader(const std::string & filename,
         const smt::SmtSolver & solver,
         MoxiEnumEncoding enum_encoding = MoxiEnumEncoding::DT);

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

  /** @return the formula that keeps a term of the sort among the values of
   *  its enumeration, if another sort encodes it that has more values; null
   *  otherwise */
  smt::Term range(const MoxiSort & sort, const smt::Term & term) const;

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
  void end_define_fun(const MoxiTerm & body, const Location & loc);

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
                          const MoxiTerm & formula,
                          const Location & loc);
  void add_subsystem(const std::string & name,
                     const std::string & system,
                     const std::vector<std::string> & args,
                     const Location & loc);
  void add_check_formula(const std::string & attribute,
                         const std::string & name,
                         const MoxiTerm & formula,
                         const Location & loc);
  void add_query(const QueryExpr & query);
  void add_queries(const std::vector<QueryExpr> & queries,
                   const Location & loc);

  // terms
  MoxiTerm make_numeral(const std::string & digits, const Location & loc);
  MoxiTerm make_decimal(const std::string & text, const Location & loc);
  MoxiTerm make_binary(const std::string & bits, const Location & loc);
  MoxiTerm make_hexadecimal(const std::string & digits, const Location & loc);
  /** @return the term an identifier stands for on its own */
  MoxiTerm make_identifier_term(const Identifier & id);
  /** @return the application of the function an identifier names */
  MoxiTerm make_application(const Identifier & id,
                            const MoxiTermVec & args,
                            const Location & loc);
  /** Binds the symbols of a let for its body. */
  void push_let(const std::vector<VarBinding> & bindings);
  /** Ends the scope of the innermost let, or of define-fun parameters. */
  void pop_scope();
  /** Binds the variables of a quantifier for its body. */
  void push_quantifier(const std::vector<SortedVar> & vars);
  /** Ends the scope of the innermost quantifier.
   *  @return the quantified formula */
  MoxiTerm pop_quantifier(bool forall,
                          const MoxiTerm & body,
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
    std::vector<MoxiSort> param_sorts;
    smt::Term body;
    MoxiSort result;
  };

  /** An uninterpreted function from declare-fun. */
  struct Function
  {
    smt::Term symbol;
    std::vector<MoxiSort> arg_sorts;
    MoxiSort result;
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

  /** A number, as the text and the sort that the solver makes it from. */
  struct Number
  {
    std::string text;
    smt::Sort sort;
  };

  /** What the command being parsed is, for the attributes it may take. */
  enum class Scope
  {
    TOP,
    SYSTEM,
    CHECK
  };

  using SortParams = std::unordered_map<std::string, MoxiSort>;

  MoxiSort resolve_sort(const SortExpr & expr,
                        const SortParams * params = nullptr);
  /** @return the one object that tells apart the sort, which involves an
   *  enumeration that another sort encodes */
  const EncodedSort * encode_sort(const smt::Sort & sort,
                                  const std::string & name,
                                  const Enumeration * enumeration,
                                  std::vector<MoxiSort> args);
  /** Refuses where the solver would tell apart values of a sort that are
   *  equal in MoXI, as the sort that encodes an enumeration has other values.
   *  @param reason what would tell them apart */
  [[noreturn]] void refuse_encoding(const Enumeration & enumeration,
                                    const std::string & reason,
                                    const Location & loc) const;
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
  const MoxiTerm * lookup_bound(const std::string & name) const;
  void bind(const std::string & name, const MoxiTerm & term);
  MoxiTerm make_next(const std::string & name, const Location & loc);

  MoxiTerm make_indexed_constant(const Identifier & id);
  MoxiTerm make_indexed_application(const Identifier & id,
                                    const MoxiTermVec & args,
                                    const Location & loc);
  MoxiTerm make_builtin(const std::string & name,
                        MoxiTermVec args,
                        const Location & loc);
  /** Applies a predefined function to arguments of which one at least
   *  involves an enumeration that another sort encodes. */
  MoxiTerm make_builtin_of_enumerations(const std::string & name,
                                        smt::PrimOp op,
                                        const MoxiTermVec & args,
                                        const Location & loc);
  MoxiTerm apply_macro(const std::string & name,
                       const Macro & macro,
                       const MoxiTermVec & args,
                       const Location & loc);
  /** @return the equality of two terms of the sort as MoXI defines it,
   *  which only the values of an enumeration take part in
   *  @param name the function that compares them, for errors */
  smt::Term make_equality(const smt::Term & a,
                          const smt::Term & b,
                          const MoxiSort & sort,
                          const std::string & name,
                          const Location & loc);
  /** @return the term, or the first value of its enumeration where it is
   *  not one of its values */
  MoxiTerm clamp(const MoxiTerm & term, const Location & loc) const;
  /** @return whether the term is a formula, i.e. has sort Bool */
  bool is_formula(const MoxiTerm & term) const;
  /** @return whether the term has the sort, or is a real with an integer
   *  value that the solver gives the sort Int */
  bool has_sort(const MoxiTerm & term, const MoxiSort & sort) const;
  /** Makes a number from its text, which numbers_ keeps. */
  smt::Term make_number(const std::string & text,
                        const smt::Sort & sort,
                        const Location & loc);
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
  MoxiEnumEncoding enum_encoding_;
  Location location_;

  // the logic
  bool logic_set_ = false;
  bool logic_has_int_ = true;
  bool logic_has_real_ = true;
  bool logic_has_bv_ = true;
  /** the sort of formulas, of kind BV for solvers that alias Bool and
   *  (_ BitVec 1), e.g. Boolector */
  smt::Sort bool_sort_;
  smt::Sort int_sort_;
  smt::Sort real_sort_;
  /** the numbers made so far, i.e. the values of numerals, decimals and
   *  their negations */
  std::unordered_map<smt::Term, Number> numbers_;

  // global declarations
  std::unordered_map<std::string, DeclaredSort> declared_sorts_;
  std::unordered_map<std::string, SortDefinition> defined_sorts_;
  std::unordered_map<std::string, MoxiSort> enum_sorts_;
  std::unordered_map<std::string, MoxiTerm> enum_values_;
  /** the enumerations that other sorts encode, and the sorts made of them;
   *  a deque keeps pointers valid */
  std::deque<Enumeration> enumerations_;
  std::deque<EncodedSort> encoded_sorts_;
  std::unordered_map<std::string, const EncodedSort *> encoded_sort_names_;
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
  std::unordered_map<std::string, std::vector<MoxiTerm>> bound_;
  std::vector<std::vector<std::string>> bound_scopes_;
  std::vector<MoxiTermVec> quantifier_params_;
};

}  // namespace moxi
}  // namespace pono
