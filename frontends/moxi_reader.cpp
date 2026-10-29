/*********************                                                        */
/*! \file moxi_reader.cpp
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
**/

#include "frontends/moxi_reader.h"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cctype>
#include <limits>
#include <ostream>
#include <sstream>

#include "utils/exceptions.h"
#include "utils/logger.h"
#include "utils/str_util.h"

using namespace smt;
using namespace std;

namespace pono {
namespace moxi {

// pono's own to_string, for its options and results, would hide it.
using std::to_string;

namespace {

/** A predefined function that maps to an operator of the solver. It takes
 *  as many arguments as smt-switch's get_arity says the operator takes. */
struct Builtin
{
  PrimOp op;
  /** a bound on the number of arguments, for a function that takes fewer
   *  than the operator */
  size_t max_args = numeric_limits<size_t>::max();
};

const unordered_map<string, Builtin> & builtins()
{
  static const unordered_map<string, Builtin> table{
    // core
    { "not", { Not } },
    { "and", { And } },
    { "or", { Or } },
    { "xor", { Xor } },
    { "=>", { Implies } },
    { "=", { Equal } },
    { "distinct", { Distinct } },
    // MoXI adds a binary disequality
    { "!=", { Distinct, 2 } },
    { "ite", { Ite } },
    // integers and reals
    { "+", { Plus } },
    { "-", { Minus } },
    { "*", { Mult } },
    { "/", { Div } },
    { "div", { IntDiv } },
    { "mod", { Mod } },
    { "abs", { Abs } },
    { "<", { Lt } },
    { "<=", { Le } },
    { ">", { Gt } },
    { ">=", { Ge } },
    { "to_real", { To_Real } },
    { "to_int", { To_Int } },
    { "is_int", { Is_Int } },
    // bit-vectors
    { "concat", { Concat } },
    { "bvnot", { BVNot } },
    { "bvneg", { BVNeg } },
    { "bvand", { BVAnd } },
    { "bvor", { BVOr } },
    { "bvxor", { BVXor } },
    { "bvadd", { BVAdd } },
    { "bvmul", { BVMul } },
    { "bvnand", { BVNand } },
    { "bvnor", { BVNor } },
    { "bvxnor", { BVXnor } },
    { "bvsub", { BVSub } },
    { "bvudiv", { BVUdiv } },
    { "bvsdiv", { BVSdiv } },
    { "bvurem", { BVUrem } },
    { "bvsrem", { BVSrem } },
    { "bvsmod", { BVSmod } },
    { "bvshl", { BVShl } },
    { "bvlshr", { BVLshr } },
    { "bvashr", { BVAshr } },
    { "bvcomp", { BVComp } },
    { "bvult", { BVUlt } },
    { "bvule", { BVUle } },
    { "bvugt", { BVUgt } },
    { "bvuge", { BVUge } },
    { "bvslt", { BVSlt } },
    { "bvsle", { BVSle } },
    { "bvsgt", { BVSgt } },
    { "bvsge", { BVSge } },
    { "ubv_to_int", { UBV_To_Int } },
    { "sbv_to_int", { SBV_To_Int } },
    // arrays
    { "select", { Select } },
    { "store", { Store } },
  };
  return table;
}

/** The operators that SMT-LIB declares as (op Real Real s), which AUFLIRA
 *  and AUFNIRA let take integers for reals. */
const unordered_set<PrimOp> real_operators{ Plus, Minus, Mult, Div,
                                            Lt,   Le,    Gt,   Ge };

/** Predefined functions that take indices, as in ((_ extract 7 0) x). */
const unordered_map<string, PrimOp> indexed_builtins{
  { "extract", Extract },
  { "zero_extend", Zero_Extend },
  { "sign_extend", Sign_Extend },
  { "rotate_left", Rotate_Left },
  { "rotate_right", Rotate_Right },
  { "repeat", Repeat },
  { "int_to_bv", Int_To_BV },
  // not an operator of the solver, see make_indexed_application
  { "divisible", NUM_OPS_AND_NULL },
};

const unordered_set<string> builtin_sorts{
  "Bool", "Int", "Real", "BitVec", "Array"
};

bool is_predefined_function(const string & name)
{
  return name == "true" || name == "false" || builtins().count(name)
         || indexed_builtins.count(name);
}

/** @return whether the name is bv followed by a numeral, as in (_ bv5 8) */
bool is_bv_literal_name(const string & name)
{
  return name.size() > 2 && name.compare(0, 2, "bv") == 0
         && all_of(name.begin() + 2, name.end(), [](char c) {
              return isdigit(static_cast<unsigned char>(c));
            });
}

/** Recognizes a logic by the SMT-LIB naming scheme: an optional QF_, the
 *  letters of other theories, and then those of the arithmetic, if any.
 *  @return whether the name is a logic
 */
bool parse_logic(const string & logic,
                 bool & has_int,
                 bool & has_real,
                 bool & has_bv)
{
  if (logic == "ALL") {
    has_int = has_real = has_bv = true;
    return true;
  }
  string rest = logic;
  if (syntax_analysis::StrStartsWith(rest, "QF_")) {
    rest = rest.substr(3);
  }
  bool has_theory = false;
  bool bv = false;
  for (const string theory : { "AX", "A", "UF", "BV", "FP", "DT", "S" }) {
    if (syntax_analysis::StrStartsWith(rest, theory)) {
      rest = rest.substr(theory.size());
      has_theory = true;
      bv |= theory == "BV";
    }
  }
  static const unordered_map<string, pair<bool, bool>> arithmetic{
    { "", { false, false } },   { "IDL", { true, false } },
    { "RDL", { false, true } }, { "LIA", { true, false } },
    { "LRA", { false, true } }, { "NIA", { true, false } },
    { "NRA", { false, true } }, { "LIRA", { true, true } },
    { "NIRA", { true, true } },
  };
  auto it = arithmetic.find(rest);
  if (it == arithmetic.end() || (rest.empty() && !has_theory)) {
    return false;
  }
  has_int = it->second.first;
  has_real = it->second.second;
  has_bv = bv;
  return true;
}

/** @return a name for a placeholder symbol. A solver refuses a second symbol
 *  of a name, and readers may share one, so the count is shared too. */
string placeholder_name()
{
  static atomic<uint64_t> count{ 0 };
  return generated_name("moxi_" + to_string(count++));
}

string describe_sorts(const TermVec & args)
{
  string result;
  for (const Term & arg : args) {
    result += (result.empty() ? "" : ", ") + arg->get_sort()->to_string();
  }
  return result;
}

string describe_sorts(const MoxiTermVec & args)
{
  string result;
  for (const MoxiTerm & arg : args) {
    result += (result.empty() ? "" : ", ") + arg.sort().to_string();
  }
  return result;
}

/** @return the sort of MoXI that is the solver's sort */
MoxiSort plain_sort(const Sort & sort) { return { sort, nullptr }; }

/** @return the term, whose sort of MoXI is its sort */
MoxiTerm plain(const Term & term) { return { term, nullptr }; }

/** @return the term with the sort of MoXI */
MoxiTerm typed(const Term & term, const MoxiSort & sort)
{
  return { term, sort.encoded };
}

MoxiSort element_sort(const MoxiSort & array)
{
  return array.encoded ? array.encoded->args[1]
                       : plain_sort(array.sort->get_elemsort());
}

/** @return the enumeration the sort is if a sort with more values encodes
 *  it, or null */
const Enumeration * inexact_enumeration(const MoxiSort & sort)
{
  const Enumeration * enumeration =
      sort.encoded ? sort.encoded->enumeration : nullptr;
  return enumeration && !enumeration->exact ? enumeration : nullptr;
}

/** @return an enumeration that a sort with more values encodes, by which the
 *  solver's equality on the sort tells apart values that are equal in MoXI,
 *  or null. An array holds elements at the other values of the encoding sort
 *  if it is indexed by the enumeration, and may hold these values if they
 *  are its elements. The other sorts have the same values in MoXI and in the
 *  solver: a term of an enumeration takes only its values, and those of a
 *  declared sort are opaque. */
const Enumeration * inexact_equality(const MoxiSort & sort)
{
  if (!sort.encoded || sort.sort->get_sort_kind() != ARRAY) {
    return nullptr;
  }
  for (const MoxiSort & arg : sort.encoded->args) {
    if (const Enumeration * enumeration = inexact_enumeration(arg)) {
      return enumeration;
    }
  }
  for (const MoxiSort & arg : sort.encoded->args) {
    if (const Enumeration * enumeration = inexact_equality(arg)) {
      return enumeration;
    }
  }
  return nullptr;
}

}  // namespace

string MoxiSort::to_string() const
{
  return encoded ? encoded->name : sort->to_string();
}

MoxiSort MoxiTerm::sort() const
{
  // the encoding sort, as the solver may report another for a value, e.g.
  // MathSAT reports Int for the real 1
  return encoded ? MoxiSort{ encoded->sort, encoded }
                 : MoxiSort{ term->get_sort(), nullptr };
}

void Location::columns(int count)
{
  const int64_t column = static_cast<int64_t>(end.column) + count;
  end.column = static_cast<uint32_t>(max<int64_t>(1, column));
}

void Location::lines(int count)
{
  if (count > 0) {
    end.line += static_cast<uint32_t>(count);
    end.column = 1;
  }
}

ostream & operator<<(ostream & os, const Location & loc)
{
  return os << loc.begin.line << ":" << loc.begin.column;
}

const Variable & System::variable(size_t i) const
{
  if (i < inputs.size()) {
    return inputs[i];
  }
  i -= inputs.size();
  if (i < outputs.size()) {
    return outputs[i];
  }
  return locals.at(i - outputs.size());
}

Reader::Reader(const string & filename,
               const SmtSolver & solver,
               MoxiEnumEncoding enum_encoding)
    : filename_(filename),
      solver_(solver),
      enum_encoding_(enum_encoding),
      bool_sort_(solver->make_sort(BOOL))
{
  parse();
}

Term Reader::range(const MoxiSort & sort, const Term & term) const
{
  const Enumeration * enumeration = inexact_enumeration(sort);
  if (!enumeration) {
    return nullptr;
  }
  const TermVec & values = enumeration->values;
  switch (enumeration->sort->get_sort_kind()) {
    case BV: return solver_->make_term(BVUle, term, values.back());
    case INT:
      return solver_->make_term(And,
                                solver_->make_term(Le, values.front(), term),
                                solver_->make_term(Le, term, values.back()));
    default: {
      // Reals lie between the values, too.
      assert(enumeration->sort->get_sort_kind() == REAL);
      TermVec equalities;
      for (const Term & value : values) {
        equalities.push_back(solver_->make_term(Equal, term, value));
      }
      return equalities.size() == 1 ? equalities[0]
                                    : solver_->make_term(Or, equalities);
    }
  }
}

void Reader::error(const Location & loc, const string & message) const
{
  ostringstream text;
  text << filename_ << ":" << loc << ": " << message;
  throw PonoException(text.str());
}

// commands

void Reader::set_logic(const string & logic, const Location & loc)
{
  if (logic_set_) {
    error(loc, "the logic is already set");
  }
  if (!parse_logic(logic, logic_has_int_, logic_has_real_, logic_has_bv_)) {
    error(loc, "unknown logic " + logic);
  }
  logic_set_ = true;
}

void Reader::declare_sort(const string & name,
                          const string & arity,
                          const Location & loc)
{
  check_new_sort_name(name, loc);
  const uint64_t n = parse_index(arity, loc);
  try {
    declared_sorts_[name] = { solver_->make_sort(name, n), n };
  }
  catch (const SmtException & e) {
    error(loc, "cannot declare sort " + name + ": " + e.what());
  }
}

void Reader::define_sort(const string & name,
                         const vector<string> & params,
                         const SortExpr & body,
                         const Location & loc)
{
  check_new_sort_name(name, loc);
  // Resolve the body once to report its mistakes here, with the parameters
  // standing for an arbitrary sort.
  SortParams dummies;
  for (const string & param : params) {
    if (!dummies.emplace(param, plain_sort(bool_sort_)).second) {
      error(loc, "sort parameter " + param + " is repeated");
    }
  }
  resolve_sort(body, &dummies);
  defined_sorts_[name] = { params, body };
}

void Reader::declare_enum_sort(const string & name,
                               const vector<string> & values,
                               const Location & loc)
{
  check_new_sort_name(name, loc);
  unordered_set<string> seen;
  for (const string & value : values) {
    check_new_global_name(value, loc);
    if (!seen.insert(value).second) {
      error(loc, "value " + value + " of enumeration " + name + " is repeated");
    }
  }
  const SolverEnum se = solver_->get_solver_enum();
  const bool has_datatypes = solver_has_attribute(se, THEORY_DATATYPE);
  MoxiEnumEncoding encoding = enum_encoding_;
  if (encoding == MoxiEnumEncoding::AUTO) {
    // Datatypes encode an enumeration exactly. Otherwise a sort of the logic
    // does, or bit-vectors, which every solver has, if the logic has none.
    if (has_datatypes) {
      encoding = MoxiEnumEncoding::DT;
    } else if (logic_has_bv_ || (!logic_has_int_ && !logic_has_real_)) {
      encoding = MoxiEnumEncoding::BV;
    } else {
      encoding =
          logic_has_int_ ? MoxiEnumEncoding::INT : MoxiEnumEncoding::REAL;
    }
  }
  if (encoding == MoxiEnumEncoding::DT && !has_datatypes) {
    error(loc,
          "solver " + to_string(se)
              + " does not support datatypes, which enumeration sorts are;"
                " choose one that does, e.g. with --smt-solver cvc5, or"
                " another encoding with --moxi-enum-as");
  }
  try {
    if (encoding == MoxiEnumEncoding::DT) {
      // An enumeration is a datatype whose constructors take no arguments.
      DatatypeDecl decl = solver_->make_datatype_decl(name);
      for (const string & value : values) {
        solver_->add_constructor(
            decl, solver_->make_datatype_constructor_decl(value));
      }
      const Sort sort = solver_->make_sort(decl);
      enum_sorts_[name] = plain_sort(sort);
      for (const string & value : values) {
        enum_values_[value] = plain(solver_->make_term(
            Apply_Constructor, solver_->get_constructor(sort, value)));
      }
      logger.log(1, "MoXI enumeration {} is a datatype", name);
      return;
    }
    // Another sort encodes the enumeration, from which its name keeps it
    // apart.
    Enumeration & enumeration = enumerations_.emplace_back();
    enumeration.name = name;
    if (encoding == MoxiEnumEncoding::BV) {
      // as few bits as tell the values apart, but at least one
      uint64_t width = 1;
      while (width < 64 && (uint64_t{ 1 } << width) < values.size()) {
        ++width;
      }
      enumeration.sort = solver_->make_sort(BV, width);
      enumeration.exact =
          width < 64 && (uint64_t{ 1 } << width) == values.size();
    } else {
      enumeration.sort =
          encoding == MoxiEnumEncoding::INT ? int_sort(loc) : real_sort(loc);
    }
    for (size_t i = 0; i < values.size(); ++i) {
      enumeration.values.push_back(
          solver_->make_term(to_string(i), enumeration.sort));
    }
    const MoxiSort sort{
      enumeration.sort, encode_sort(enumeration.sort, name, &enumeration, {})
    };
    enum_sorts_[name] = sort;
    for (size_t i = 0; i < values.size(); ++i) {
      enum_values_[values[i]] = typed(enumeration.values[i], sort);
    }
    logger.log(1,
               "MoXI enumeration {} is encoded as {}",
               name,
               enumeration.sort->to_string());
  }
  catch (const SmtException & e) {
    error(loc, "cannot declare enumeration " + name + ": " + e.what());
  }
}

void Reader::declare_fun(const string & name,
                         const vector<SortExpr> & args,
                         const SortExpr & result,
                         const Location & loc)
{
  check_new_global_name(name, loc);
  const MoxiSort result_sort = resolve_sort(result);
  if (args.empty()) {
    constants_.push_back(
        { name, result_sort, make_placeholder(result_sort.sort), loc });
    constant_names_[name] = &constants_.back();
    return;
  }
  Function function;
  SortVec sorts;
  for (const SortExpr & arg : args) {
    function.arg_sorts.push_back(resolve_sort(arg));
    sorts.push_back(function.arg_sorts.back().sort);
    // A function gives equal results for equal arguments, but the solver
    // would tell apart some that are equal in MoXI.
    if (const Enumeration * enumeration =
            inexact_equality(function.arg_sorts.back())) {
      refuse_encoding(*enumeration,
                      "function " + name
                          + " would tell apart equal arguments of sort "
                          + function.arg_sorts.back().to_string(),
                      arg.loc);
    }
  }
  sorts.push_back(result_sort.sort);
  function.result = result_sort;
  try {
    function.symbol =
        solver_->make_symbol(name, solver_->make_sort(FUNCTION, sorts));
  }
  catch (const SmtException & e) {
    error(loc, "cannot declare function " + name + ": " + e.what());
  }
  functions_[name] = function;
  function_names_.insert(name);
}

void Reader::begin_define_fun(const string & name,
                              const vector<SortedVar> & params,
                              const SortExpr & result,
                              const Location & loc)
{
  check_new_global_name(name, loc);
  macro_name_ = name;
  macro_ = Macro();
  bound_scopes_.emplace_back();
  unordered_set<string> seen;
  for (const SortedVar & param : params) {
    if (!seen.insert(param.name).second) {
      error(param.loc, "parameter " + param.name + " is repeated");
    }
    const MoxiSort sort = resolve_sort(*param.sort);
    macro_.params.push_back(make_placeholder(sort.sort));
    macro_.param_sorts.push_back(sort);
    bind(param.name, typed(macro_.params.back(), sort));
  }
  macro_.result = resolve_sort(result);
  next_allowed_ = false;
  formula_attribute_ = "define-fun";
}

void Reader::end_define_fun(const MoxiTerm & body, const Location & loc)
{
  pop_scope();
  macro_.body = body.term;
  if (!has_sort(body, macro_.result)) {
    error(loc,
          "the definition of " + macro_name_ + " has sort "
              + body.sort().to_string() + ", but " + macro_.result.to_string()
              + " is declared");
  }
  macros_[macro_name_] = std::move(macro_);
}

void Reader::begin_system(const string & name, const Location & loc)
{
  scope_ = Scope::SYSTEM;
  system_ = make_shared<System>();
  system_->name = name;
  system_->loc = loc;
  declared_attributes_.clear();
  command_names_.clear();
  variables_in_scope_.clear();
}

void Reader::end_system()
{
  // A later definition of the same name replaces this one, while the
  // systems defined in between keep instantiating this one.
  systems_[system_->name] = system_;
  system_.reset();
  variables_in_scope_.clear();
  scope_ = Scope::TOP;
}

void Reader::begin_check(const string & system, const Location & loc)
{
  auto it = systems_.find(system);
  if (it == systems_.end()) {
    error(loc, "unknown system " + system);
  }
  scope_ = Scope::CHECK;
  check_ = Check();
  check_.system = it->second;
  check_.loc = loc;
  declared_attributes_.clear();
  command_names_.clear();
  variables_in_scope_.clear();
  formula_positions_.clear();
  query_exprs_.clear();
}

void Reader::end_check()
{
  unordered_set<string> query_names;
  for (const QueryExpr & expr : query_exprs_) {
    if (!query_names.insert(expr.name).second) {
      error(expr.loc, "query " + expr.name + " is defined more than once");
    }
    if (expr.formulas.empty()) {
      error(expr.loc, "query " + expr.name + " lists no formulas");
    }
    Query query{ expr.name, {}, expr.loc };
    size_t currents = 0;
    for (const string & name : expr.formulas) {
      auto it = formula_positions_.find(name);
      if (it == formula_positions_.end()) {
        error(expr.loc,
              "query " + expr.name + " refers to " + name
                  + ", which is not a formula of this check-system");
      }
      if (find(query.formulas.begin(), query.formulas.end(), it->second)
          != query.formulas.end()) {
        continue;
      }
      query.formulas.push_back(it->second);
      if (check_.formulas[it->second].kind == FormulaKind::CURRENT) {
        ++currents;
      }
    }
    if (currents > 1) {
      error(expr.loc,
            "query " + expr.name + " has more than one :current formula");
    }
    check_.queries.push_back(std::move(query));
  }
  checks_.push_back(std::move(check_));
  check_ = Check();
  variables_in_scope_.clear();
  scope_ = Scope::TOP;
}

vector<Variable> & Reader::variables(const string & attribute,
                                     const Location & loc)
{
  if (scope_ == Scope::SYSTEM) {
    if (attribute == ":input") return system_->inputs;
    if (attribute == ":output") return system_->outputs;
    if (attribute == ":local") return system_->locals;
  } else if (scope_ == Scope::CHECK) {
    if (attribute == ":input") return check_.inputs;
    if (attribute == ":output") return check_.outputs;
    if (attribute == ":local") return check_.locals;
  }
  error(loc, "unexpected attribute " + attribute);
}

void Reader::declare_variables(const string & attribute,
                               const vector<SortedVar> & vars,
                               const Location & loc)
{
  if (!declared_attributes_.insert(attribute).second) {
    error(loc, "attribute " + attribute + " appears more than once");
  }
  vector<Variable> & group = variables(attribute, loc);
  for (const SortedVar & var : vars) {
    if (is_generated_name(var.name)) {
      error(var.loc, "the name " + var.name + " is reserved for pono");
    }
    if (!command_names_.insert(var.name).second) {
      error(var.loc, "variable " + var.name + " is declared more than once");
    }
    const MoxiSort sort = resolve_sort(*var.sort);
    group.push_back({ var.name,
                      sort,
                      make_placeholder(sort.sort),
                      make_placeholder(sort.sort),
                      var.loc });
  }
}

void Reader::check_signature(const string & attribute,
                             vector<Variable> & declared,
                             const vector<Variable> & expected,
                             const Location & loc)
{
  const System & system = *check_.system;
  if (!declared_attributes_.count(attribute)) {
    // Leaving out the attribute keeps the names the system gave them.
    for (const Variable & var : expected) {
      if (!command_names_.insert(var.name).second) {
        error(loc,
              "variable " + var.name + " of system " + system.name
                  + " clashes with a variable of this check-system");
      }
      declared.push_back({ var.name,
                           var.sort,
                           make_placeholder(var.sort.sort),
                           make_placeholder(var.sort.sort),
                           loc });
    }
    return;
  }
  if (declared.size() != expected.size()) {
    error(loc,
          "the " + attribute + " attribute declares "
              + to_string(declared.size()) + " variables, but system "
              + system.name + " has " + to_string(expected.size()));
  }
  for (size_t i = 0; i < declared.size(); ++i) {
    if (declared[i].sort != expected[i].sort) {
      error(declared[i].loc,
            "variable " + declared[i].name + " has sort "
                + declared[i].sort.to_string() + ", but the variable "
                + expected[i].name + " of system " + system.name
                + " it renames has sort " + expected[i].sort.to_string());
    }
  }
}

void Reader::enter_variables(const vector<Variable> & vars)
{
  for (const Variable & var : vars) {
    variables_in_scope_[var.name] = &var;
  }
}

void Reader::end_variables()
{
  variables_in_scope_.clear();
  variable_positions_.clear();
  if (scope_ == Scope::CHECK) {
    const System & system = *check_.system;
    check_signature(":input", check_.inputs, system.inputs, check_.loc);
    check_signature(":output", check_.outputs, system.outputs, check_.loc);
    check_signature(":local", check_.locals, system.locals, check_.loc);
    enter_variables(check_.inputs);
    enter_variables(check_.outputs);
    enter_variables(check_.locals);
    return;
  }
  assert(scope_ == Scope::SYSTEM);
  enter_variables(system_->inputs);
  enter_variables(system_->outputs);
  enter_variables(system_->locals);
  for (size_t i = 0; i < system_->num_variables(); ++i) {
    variable_positions_[system_->variable(i).name] = i;
  }
}

void Reader::begin_formula(const string & attribute)
{
  next_allowed_ = attribute == ":trans" || attribute == ":assumption"
                  || attribute == ":fairness" || attribute == ":reachable";
  formula_attribute_ = attribute;
}

void Reader::set_system_formula(const string & attribute,
                                const MoxiTerm & formula,
                                const Location & loc)
{
  next_allowed_ = false;
  if (!declared_attributes_.insert(attribute).second) {
    error(loc, "attribute " + attribute + " appears more than once");
  }
  if (!is_formula(formula)) {
    error(loc,
          "the formula of " + attribute + " has sort "
              + formula.sort().to_string() + " instead of Bool");
  }
  if (attribute == ":init") {
    system_->init = formula.term;
  } else if (attribute == ":trans") {
    system_->trans = formula.term;
  } else {
    assert(attribute == ":inv");
    system_->inv = formula.term;
  }
}

void Reader::add_subsystem(const string & name,
                           const string & system,
                           const vector<string> & args,
                           const Location & loc)
{
  // The instance names the variables it has, see MoxiEncoder.
  if (is_generated_name(name)) {
    error(loc, "the name " + name + " is reserved for pono");
  }
  if (!command_names_.insert(name).second) {
    error(loc,
          "subsystem " + name + " has the name of another variable or"
              + " subsystem of " + system_->name);
  }
  auto it = systems_.find(system);
  if (it == systems_.end()) {
    error(loc, "unknown system " + system);
  }
  const System & target = *it->second;
  const size_t expected = target.inputs.size() + target.outputs.size();
  if (args.size() != expected) {
    error(loc,
          "subsystem " + name + " passes " + to_string(args.size())
              + " variables to system " + system + ", which has "
              + to_string(expected) + " inputs and outputs");
  }
  Subsystem subsystem{ name, it->second, {}, loc };
  for (size_t i = 0; i < args.size(); ++i) {
    auto pos = variable_positions_.find(args[i]);
    if (pos == variable_positions_.end()) {
      error(loc,
            "subsystem " + name + " passes " + args[i]
                + ", which is not a variable of system " + system_->name);
    }
    const Variable & actual = system_->variable(pos->second);
    const Variable & formal = target.variable(i);
    if (actual.sort != formal.sort) {
      error(loc,
            "subsystem " + name + " passes " + actual.name + " of sort "
                + actual.sort.to_string() + " for " + formal.name + " of sort "
                + formal.sort.to_string());
    }
    subsystem.args.push_back(pos->second);
  }
  system_->subsystems.push_back(std::move(subsystem));
}

void Reader::add_check_formula(const string & attribute,
                               const string & name,
                               const MoxiTerm & formula,
                               const Location & loc)
{
  next_allowed_ = false;
  if (!is_formula(formula)) {
    error(loc,
          "the formula " + name + " has sort " + formula.sort().to_string()
              + " instead of Bool");
  }
  FormulaKind kind = FormulaKind::CURRENT;
  if (attribute == ":assumption") {
    kind = FormulaKind::ASSUMPTION;
  } else if (attribute == ":fairness") {
    kind = FormulaKind::FAIRNESS;
  } else if (attribute == ":reachable") {
    kind = FormulaKind::REACHABLE;
  } else {
    assert(attribute == ":current");
  }
  if (!formula_positions_.emplace(name, check_.formulas.size()).second) {
    error(loc, "another formula of this check-system is named " + name);
  }
  check_.formulas.push_back({ name, kind, formula.term, loc });
}

void Reader::add_query(const QueryExpr & query)
{
  // Queries may name formulas that come later, so they are resolved at the
  // end of the command.
  query_exprs_.push_back(query);
}

void Reader::add_queries(const vector<QueryExpr> & queries,
                         const Location & loc)
{
  check_.queries_attributes.push_back(loc);
  for (const QueryExpr & query : queries) {
    add_query(query);
  }
}

// sorts

Sort Reader::int_sort(const Location & loc)
{
  if (!int_sort_) {
    const SolverEnum se = solver_->get_solver_enum();
    if (!solver_has_attribute(se, THEORY_INT)) {
      error(loc,
            "solver " + to_string(se)
                + " does not support integers; choose one that does,"
                  " e.g. with --smt-solver cvc5");
    }
    int_sort_ = solver_->make_sort(INT);
  }
  return int_sort_;
}

Sort Reader::real_sort(const Location & loc)
{
  if (!real_sort_) {
    const SolverEnum se = solver_->get_solver_enum();
    if (!solver_has_attribute(se, THEORY_REAL)) {
      error(loc,
            "solver " + to_string(se)
                + " does not support reals; choose one that does,"
                  " e.g. with --smt-solver cvc5");
    }
    real_sort_ = solver_->make_sort(REAL);
  }
  return real_sort_;
}

Sort Reader::numeral_sort(const Location & loc)
{
  // In SMT-LIB, numerals are reals in the logics without integers.
  return logic_has_real_ && !logic_has_int_ ? real_sort(loc) : int_sort(loc);
}

uint64_t Reader::parse_index(const string & index, const Location & loc) const
{
  if (index.empty() || index.size() > 19
      || !all_of(index.begin(), index.end(), [](char c) {
           return isdigit(static_cast<unsigned char>(c));
         })) {
    error(loc, "expected a numeral that fits 64 bits, not " + index);
  }
  return stoull(index);
}

MoxiSort Reader::resolve_sort(const SortExpr & expr, const SortParams * params)
{
  const string & name = expr.name;
  const bool plain = expr.indices.empty() && expr.args.empty();
  if (params && plain) {
    auto it = params->find(name);
    if (it != params->end()) {
      return it->second;
    }
  }
  try {
    if (name == "Bool" && plain) {
      return plain_sort(bool_sort_);
    } else if (name == "Int" && plain) {
      return plain_sort(int_sort(expr.loc));
    } else if (name == "Real" && plain) {
      return plain_sort(real_sort(expr.loc));
    } else if (name == "BitVec" && expr.indices.size() == 1
               && expr.args.empty()) {
      const uint64_t width = parse_index(expr.indices[0], expr.loc);
      if (!width) {
        error(expr.loc, "bit-vectors need a positive width");
      }
      return plain_sort(solver_->make_sort(BV, width));
    } else if (name == "Array" && expr.indices.empty()
               && expr.args.size() == 2) {
      const MoxiSort index = resolve_sort(expr.args[0], params);
      const MoxiSort element = resolve_sort(expr.args[1], params);
      // Selecting at equal indices gives equal elements, but the solver would
      // tell apart some that are equal in MoXI.
      if (const Enumeration * enumeration = inexact_equality(index)) {
        refuse_encoding(*enumeration,
                        "an array would tell apart equal indices of sort "
                            + index.to_string(),
                        expr.args[0].loc);
      }
      MoxiSort array =
          plain_sort(solver_->make_sort(ARRAY, index.sort, element.sort));
      if (index.involves_enumeration() || element.involves_enumeration()) {
        array.encoded = encode_sort(
            array.sort,
            "(Array " + index.to_string() + " " + element.to_string() + ")",
            nullptr,
            { index, element });
      }
      return array;
    }

    auto enumeration = enum_sorts_.find(name);
    if (enumeration != enum_sorts_.end() && plain) {
      return enumeration->second;
    }

    auto defined = defined_sorts_.find(name);
    if (defined != defined_sorts_.end() && expr.indices.empty()
        && expr.args.size() == defined->second.params.size()) {
      SortParams args;
      for (size_t i = 0; i < expr.args.size(); ++i) {
        args[defined->second.params[i]] = resolve_sort(expr.args[i], params);
      }
      return resolve_sort(defined->second.body, &args);
    }

    auto declared = declared_sorts_.find(name);
    if (declared != declared_sorts_.end() && expr.indices.empty()
        && expr.args.size() == declared->second.arity) {
      if (!declared->second.arity) {
        return plain_sort(declared->second.sort);
      }
      vector<MoxiSort> args;
      SortVec sorts;
      for (const SortExpr & arg : expr.args) {
        args.push_back(resolve_sort(arg, params));
        sorts.push_back(args.back().sort);
      }
      MoxiSort sort =
          plain_sort(solver_->make_sort(declared->second.sort, sorts));
      if (any_of(args.begin(), args.end(), [](const MoxiSort & arg) {
            return arg.involves_enumeration();
          })) {
        string text = "(" + name;
        for (const MoxiSort & arg : args) {
          text += " " + arg.to_string();
        }
        sort.encoded = encode_sort(sort.sort, text + ")", nullptr, args);
      }
      return sort;
    }
  }
  catch (const SmtException & e) {
    error(expr.loc, "cannot make sort " + name + ": " + e.what());
  }

  const bool known = builtin_sorts.count(name) || enum_sorts_.count(name)
                     || defined_sorts_.count(name)
                     || declared_sorts_.count(name);
  error(expr.loc,
        known ? "wrong indices or arguments for sort " + name
              : "unknown sort " + name);
}

const EncodedSort * Reader::encode_sort(const Sort & sort,
                                        const string & name,
                                        const Enumeration * enumeration,
                                        vector<MoxiSort> args)
{
  // How MoXI writes the sort tells it apart, as enumerations have names of
  // their own.
  auto it = encoded_sort_names_.find(name);
  if (it != encoded_sort_names_.end()) {
    return it->second;
  }
  encoded_sorts_.push_back({ sort, name, enumeration, std::move(args) });
  return encoded_sort_names_[name] = &encoded_sorts_.back();
}

void Reader::refuse_encoding(const Enumeration & enumeration,
                             const string & reason,
                             const Location & loc) const
{
  error(loc,
        "cannot encode enumeration " + enumeration.name + " as "
            + enumeration.sort->to_string() + " here, as " + reason
            + " by values outside " + enumeration.name
            + "; datatypes (--moxi-enum-as dt) encode it exactly");
}

void Reader::check_new_sort_name(const string & name,
                                 const Location & loc) const
{
  if (builtin_sorts.count(name) || enum_sorts_.count(name)
      || defined_sorts_.count(name) || declared_sorts_.count(name)) {
    error(loc, "sort " + name + " is already defined");
  }
}

void Reader::check_new_global_name(const string & name,
                                   const Location & loc) const
{
  if (is_generated_name(name)) {
    error(loc, "the name " + name + " is reserved for pono");
  }
  if (is_predefined_function(name)) {
    error(loc, name + " is a predefined symbol");
  }
  if (enum_values_.count(name) || constant_names_.count(name)
      || functions_.count(name) || macros_.count(name)) {
    error(loc, "symbol " + name + " is already declared");
  }
}

Term Reader::make_placeholder(const Sort & sort)
{
  return solver_->make_symbol(placeholder_name(), sort);
}

// terms

const Variable * Reader::lookup_variable(const string & name) const
{
  auto it = variables_in_scope_.find(name);
  return it == variables_in_scope_.end() ? nullptr : it->second;
}

const MoxiTerm * Reader::lookup_bound(const string & name) const
{
  auto it = bound_.find(name);
  return it == bound_.end() ? nullptr : &it->second.back();
}

void Reader::bind(const string & name, const MoxiTerm & term)
{
  bound_[name].push_back(term);
  bound_scopes_.back().push_back(name);
}

void Reader::push_let(const vector<VarBinding> & bindings)
{
  // The bindings of a let are parallel: each term was built without the
  // others in scope, so they are all bound only now.
  bound_scopes_.emplace_back();
  unordered_set<string> seen;
  for (const VarBinding & binding : bindings) {
    if (!seen.insert(binding.name).second) {
      error(binding.loc, "let binds " + binding.name + " more than once");
    }
    bind(binding.name, binding.term);
  }
}

void Reader::pop_scope()
{
  for (const string & name : bound_scopes_.back()) {
    // bind() put the name in, so this does not insert it
    vector<MoxiTerm> & terms = bound_[name];
    assert(!terms.empty());
    if (terms.size() == 1) {
      bound_.erase(name);
    } else {
      terms.pop_back();
    }
  }
  bound_scopes_.pop_back();
}

void Reader::push_quantifier(const vector<SortedVar> & vars)
{
  bound_scopes_.emplace_back();
  MoxiTermVec params;
  unordered_set<string> seen;
  for (const SortedVar & var : vars) {
    if (!seen.insert(var.name).second) {
      error(var.loc, "quantifier binds " + var.name + " more than once");
    }
    const MoxiSort sort = resolve_sort(*var.sort);
    // The solver's variable has the encoding sort, and keeps the sort of
    // MoXI for the guard that pop_quantifier adds.
    params.push_back(
        typed(solver_->make_param(placeholder_name(), sort.sort), sort));
    bind(var.name, params.back());
  }
  quantifier_params_.push_back(std::move(params));
}

MoxiTerm Reader::pop_quantifier(bool forall,
                                const MoxiTerm & body,
                                const Location & loc)
{
  const MoxiTermVec params = std::move(quantifier_params_.back());
  quantifier_params_.pop_back();
  pop_scope();
  // The solver refuses a body that is not a formula, unless an enumeration
  // passes for Bool with it, as two values in (_ BitVec 1) do with Boolector.
  if (body.encoded) {
    error(loc,
          "the body of the quantifier has sort " + body.sort().to_string()
              + " instead of Bool");
  }
  // In MoXI, a variable of an enumeration ranges over its values. Bound in
  // the solver, it ranges over the whole sort that encodes them, e.g. all
  // integers, so the quantifier is restricted to the values that stand for
  // those of the enumeration:
  //   (forall ((c Color)) phi) becomes (forall ((c Int)) (=> in_range phi))
  //   (exists ((c Color)) phi) becomes (exists ((c Int)) (and in_range phi))
  // where in_range is (and (<= 0 c) (<= c 2)) for three values, see range.
  // As in_range holds for exactly the values that stand for those of the
  // enumeration, the formula keeps its meaning. The guard conjoins the
  // ranges of all variables that need one; variables of other sorts, and of
  // enumerations that datatypes or 2^k values of bit-vectors encode, need
  // none.
  TermVec args;
  TermVec ranges;
  for (const MoxiTerm & param : params) {
    args.push_back(param.term);
    if (const Term in_range = range(param.sort(), param.term)) {
      ranges.push_back(in_range);
    }
  }
  Term formula = body.term;
  if (!ranges.empty()) {
    const Term guard =
        ranges.size() == 1 ? ranges[0] : make_term(And, ranges, loc);
    formula = make_term(forall ? Implies : And, { guard, formula }, loc);
  }
  args.push_back(formula);
  return plain(make_term(forall ? Forall : Exists, args, loc));
}

MoxiTerm Reader::make_numeral(const string & digits, const Location & loc)
{
  return plain(make_number(digits, numeral_sort(loc), loc));
}

MoxiTerm Reader::make_decimal(const string & text, const Location & loc)
{
  // Solvers take a decimal as a fraction more readily than as is.
  const size_t point = text.find('.');
  string denominator = "1" + string(text.size() - point - 1, '0');
  string numerator = text.substr(0, point) + text.substr(point + 1);
  numerator.erase(0,
                  min(numerator.find_first_not_of('0'), numerator.size() - 1));
  return plain(make_number(numerator + "/" + denominator, real_sort(loc), loc));
}

MoxiTerm Reader::make_binary(const string & bits, const Location & loc)
{
  try {
    return plain(
        solver_->make_term(bits, solver_->make_sort(BV, bits.size()), 2));
  }
  catch (const SmtException & e) {
    error(loc, "cannot make bit-vector #b" + bits + ": " + e.what());
  }
}

MoxiTerm Reader::make_hexadecimal(const string & digits, const Location & loc)
{
  try {
    return plain(solver_->make_term(
        digits, solver_->make_sort(BV, 4 * digits.size()), 16));
  }
  catch (const SmtException & e) {
    error(loc, "cannot make bit-vector #x" + digits + ": " + e.what());
  }
}

MoxiTerm Reader::make_next(const string & name, const Location & loc)
{
  if (lookup_bound(name)) {
    error(loc, "the bound variable " + name + " cannot be primed");
  }
  const Variable * var = lookup_variable(name);
  if (!var) {
    if (enum_values_.count(name) || constant_names_.count(name)
        || macros_.count(name) || functions_.count(name)) {
      error(loc, "only variables of a system can be primed, unlike " + name);
    }
    error(loc, "unknown variable " + name + " (primed)");
  }
  if (!next_allowed_) {
    error(loc,
          "the next-state variable " + name + "' cannot appear in "
              + formula_attribute_);
  }
  return typed(var->next, var->sort);
}

MoxiTerm Reader::make_identifier_term(const Identifier & id)
{
  if (!id.indices.empty()) {
    return make_indexed_constant(id);
  }
  if (id.qualifier) {
    if (id.name == "const") {
      error(id.loc, "a constant array needs the value of its elements");
    }
    Identifier plain = id;
    plain.qualifier.reset();
    const MoxiSort sort = resolve_sort(*id.qualifier);
    const MoxiTerm term = make_identifier_term(plain);
    if (!has_sort(term, sort)) {
      error(id.loc, id.name + " does not have sort " + sort.to_string());
    }
    return term;
  }
  if (id.primed) {
    return make_next(id.name, id.loc);
  }

  const string & name = id.name;
  if (name == "true" || name == "false") {
    return plain(solver_->make_term(name == "true"));
  }
  if (const MoxiTerm * term = lookup_bound(name)) {
    return *term;
  }
  if (const Variable * var = lookup_variable(name)) {
    return typed(var->curr, var->sort);
  }
  auto value = enum_values_.find(name);
  if (value != enum_values_.end()) {
    return value->second;
  }
  auto constant = constant_names_.find(name);
  if (constant != constant_names_.end()) {
    return typed(constant->second->placeholder, constant->second->sort);
  }
  auto macro = macros_.find(name);
  if (macro != macros_.end()) {
    if (!macro->second.params.empty()) {
      error(id.loc, "function " + name + " needs arguments");
    }
    return typed(macro->second.body, macro->second.result);
  }
  if (functions_.count(name) || is_predefined_function(name)) {
    error(id.loc, "function " + name + " needs arguments");
  }
  // A quoted symbol may carry the prime inside, as in |x'|.
  if (name.size() > 1 && name.back() == '\''
      && lookup_variable(name.substr(0, name.size() - 1))) {
    return make_next(name.substr(0, name.size() - 1), id.loc);
  }
  error(id.loc, "unknown symbol " + name);
}

MoxiTerm Reader::make_indexed_constant(const Identifier & id)
{
  if (!is_bv_literal_name(id.name) || id.indices.size() != 1 || id.primed
      || id.qualifier) {
    error(id.loc, "unknown indexed constant " + id.name);
  }
  const uint64_t width = parse_index(id.indices[0], id.loc);
  if (!width) {
    error(id.loc, "bit-vectors need a positive width");
  }
  try {
    return plain(
        solver_->make_term(id.name.substr(2), solver_->make_sort(BV, width)));
  }
  catch (const SmtException & e) {
    error(id.loc, "cannot make bit-vector " + id.name + ": " + e.what());
  }
}

MoxiTerm Reader::make_application(const Identifier & id,
                                  const MoxiTermVec & args,
                                  const Location & loc)
{
  if (id.primed) {
    error(id.loc, "the primed symbol " + id.name + "' is not a function");
  }
  if (id.qualifier) {
    const MoxiSort sort = resolve_sort(*id.qualifier);
    if (id.name == "const" && id.indices.empty()) {
      if (sort.sort->get_sort_kind() != ARRAY) {
        error(id.loc, "(as const S) needs an array sort S");
      }
      if (args.size() != 1) {
        error(loc, "a constant array takes a single value");
      }
      const MoxiSort element = element_sort(sort);
      if (!has_sort(args[0], element)) {
        error(loc,
              "a constant array of sort " + sort.to_string()
                  + " needs a value of sort " + element.to_string() + ", not "
                  + args[0].sort().to_string());
      }
      try {
        return typed(solver_->make_term(args[0].term, sort.sort), sort);
      }
      catch (const SmtException & e) {
        error(loc, string("cannot make constant array: ") + e.what());
      }
    }
    // Otherwise the sort only disambiguates, which the arguments do already.
    Identifier plain = id;
    plain.qualifier.reset();
    const MoxiTerm term = make_application(plain, args, loc);
    if (!has_sort(term, sort)) {
      error(loc,
            "the application of " + id.name + " does not have sort "
                + sort.to_string());
    }
    return term;
  }
  if (!id.indices.empty()) {
    return make_indexed_application(id, args, loc);
  }

  // Declared functions come first, as a variable cannot be applied anyway.
  auto macro = macros_.find(id.name);
  if (macro != macros_.end()) {
    return apply_macro(id.name, macro->second, args, loc);
  }
  auto function = functions_.find(id.name);
  if (function != functions_.end()) {
    const Function & f = function->second;
    if (args.size() != f.arg_sorts.size()) {
      error(loc,
            "function " + id.name + " takes " + to_string(f.arg_sorts.size())
                + " arguments, not " + to_string(args.size()));
    }
    TermVec children{ f.symbol };
    for (size_t i = 0; i < args.size(); ++i) {
      // The solver checks the other arguments, but cannot tell an
      // enumeration from the sort that encodes it.
      if ((args[i].encoded || f.arg_sorts[i].involves_enumeration())
          && !has_sort(args[i], f.arg_sorts[i])) {
        error(loc,
              "argument " + to_string(i + 1) + " of " + id.name + " has sort "
                  + args[i].sort().to_string() + " instead of "
                  + f.arg_sorts[i].to_string());
      }
      children.push_back(args[i].term);
    }
    // The function may give any value of the encoding sort, which stands
    // for a value of the enumeration.
    return clamp(typed(make_term(Apply, children, loc), f.result), loc);
  }
  return make_builtin(id.name, args, loc);
}

MoxiTerm Reader::apply_macro(const string & name,
                             const Macro & macro,
                             const MoxiTermVec & args,
                             const Location & loc)
{
  if (args.size() != macro.params.size()) {
    error(loc,
          "function " + name + " takes " + to_string(macro.params.size())
              + " arguments, not " + to_string(args.size()));
  }
  UnorderedTermMap substitution;
  for (size_t i = 0; i < args.size(); ++i) {
    const MoxiTerm & arg = args[i];
    if (!has_sort(arg, macro.param_sorts[i])) {
      error(loc,
            "argument " + to_string(i + 1) + " of " + name + " has sort "
                + arg.sort().to_string() + " instead of "
                + macro.param_sorts[i].to_string());
    }
    substitution[macro.params[i]] = arg.term;
  }
  try {
    return typed(solver_->substitute(macro.body, substitution), macro.result);
  }
  catch (const SmtException & e) {
    error(loc, "cannot apply " + name + ": " + e.what());
  }
}

MoxiTerm Reader::make_indexed_application(const Identifier & id,
                                          const MoxiTermVec & args,
                                          const Location & loc)
{
  auto it = indexed_builtins.find(id.name);
  if (it == indexed_builtins.end()) {
    error(id.loc, "unknown indexed function " + id.name);
  }
  vector<uint64_t> indices;
  for (const string & index : id.indices) {
    indices.push_back(parse_index(index, id.loc));
  }
  const size_t num_indices = id.name == "extract" ? 2 : 1;
  if (indices.size() != num_indices) {
    error(id.loc,
          id.name + " takes " + to_string(num_indices) + " indices, not "
              + to_string(indices.size()));
  }
  if (args.size() != 1) {
    error(loc,
          id.name + " takes a single argument, not " + to_string(args.size()));
  }
  if (args[0].encoded) {
    error(loc,
          "cannot apply " + id.name + " to an argument of sort "
              + args[0].sort().to_string());
  }
  if (id.name == "divisible") {
    if (!indices[0]) {
      error(id.loc, "divisible needs a positive index");
    }
    const Term divisor =
        solver_->make_term(to_string(indices[0]), int_sort(id.loc));
    const Term zero = solver_->make_term(int64_t{ 0 }, int_sort(id.loc));
    return plain(make_term(
        Equal, { make_term(Mod, { args[0].term, divisor }, loc), zero }, loc));
  }
  const Op op = num_indices == 2 ? Op(it->second, indices[0], indices[1])
                                 : Op(it->second, indices[0]);
  return plain(make_term(op, { args[0].term }, loc));
}

MoxiTerm Reader::make_builtin(const string & name,
                              MoxiTermVec args,
                              const Location & loc)
{
  auto it = builtins().find(name);
  if (it == builtins().end()) {
    error(loc, "unknown function " + name);
  }
  const Builtin & builtin = it->second;
  // SMT-LIB declares - also as negation, which smt-switch calls Negate
  const PrimOp op =
      builtin.op == Minus && args.size() == 1 ? Negate : builtin.op;
  const pair<size_t, size_t> arity = get_arity(op);
  if (args.size() < arity.first
      || args.size() > min(arity.second, builtin.max_args)) {
    error(loc,
          name + " cannot take " + to_string(args.size()) + " argument"
              + (args.size() == 1 ? "" : "s"));
  }
  if (any_of(args.begin(), args.end(), [](const MoxiTerm & arg) {
        return arg.encoded;
      })) {
    return make_builtin_of_enumerations(name, op, args, loc);
  }
  TermVec terms;
  for (const MoxiTerm & arg : args) {
    terms.push_back(arg.term);
  }

  // In a logic with both integers and reals, an integer argument of an
  // operator that SMT-LIB declares as (op Real Real s) stands for a real
  // next to a real argument, and always for /, as AUFLIRA and AUFNIRA
  // define. So does an integer argument of =, distinct and != next to a real
  // one, as cvc5 reads them outside its strict mode. Converting each integer
  // argument gives the value of the conversions in the binary applications
  // that SMT-LIB expands an n-ary one into.
  if (logic_has_int_ && logic_has_real_
      && (real_operators.count(op) || op == Equal || op == Distinct)
      && (op == Div || any_of(terms.begin(), terms.end(), [](const Term & arg) {
            return arg->get_sort()->get_sort_kind() == REAL;
          }))) {
    for (Term & arg : terms) {
      if (arg->get_sort()->get_sort_kind() == INT) {
        arg = make_term(To_Real, { arg }, loc);
      }
    }
  }

  if (op == Negate) {
    // The negation of a number stays a value, as some solvers, e.g. cvc5,
    // take only values as the elements of constant arrays. It takes the sort
    // of the number, which a solver may not report, e.g. MathSAT makes 1.0 an
    // integer.
    auto number = numbers_.find(terms[0]);
    if (number != numbers_.end()) {
      const string & text = number->second.text;
      return plain(make_number(text[0] == '-' ? text.substr(1) : "-" + text,
                               number->second.sort,
                               loc));
    }
  }
  // smt-switch applies an operator that SMT-LIB defines as :left-assoc,
  // :right-assoc, :chainable or :pairwise to any number of arguments, as
  // SMT-LIB defines it.
  return plain(make_term(op, terms, loc));
}

MoxiTerm Reader::make_builtin_of_enumerations(const string & name,
                                              PrimOp op,
                                              const MoxiTermVec & args,
                                              const Location & loc)
{
  // An enumeration has only the functions that MoXI defines for every sort,
  // which the solver cannot tell from those of the sort that encodes it.
  const string refusal =
      "cannot apply " + name + " to arguments of sort " + describe_sorts(args);
  TermVec terms;
  for (const MoxiTerm & arg : args) {
    terms.push_back(arg.term);
  }
  switch (op) {
    case Equal:
    case Distinct: {
      const MoxiSort sort = args[0].sort();
      for (const MoxiTerm & arg : args) {
        if (!has_sort(arg, sort)) {
          error(loc, refusal);
        }
      }
      // = holds if each argument equals the next, and distinct if no two
      // arguments are equal, by the equality of MoXI.
      TermVec formulas;
      for (size_t i = 0; i + 1 < args.size(); ++i) {
        for (size_t j = i + 1; j < (op == Equal ? i + 2 : args.size()); ++j) {
          const Term equal =
              make_equality(args[i].term, args[j].term, sort, name, loc);
          formulas.push_back(op == Equal ? equal
                                         : make_term(Not, { equal }, loc));
        }
      }
      return plain(formulas.size() == 1 ? formulas[0]
                                        : make_term(And, formulas, loc));
    }
    case Ite:
      if (!is_formula(args[0]) || !has_sort(args[2], args[1].sort())) {
        error(loc, refusal);
      }
      return typed(make_term(Ite, terms, loc), args[1].sort());
    case Select:
    case Store: {
      const MoxiSort array = args[0].sort();
      if (array.sort->get_sort_kind() != ARRAY || !array.encoded
          || !has_sort(args[1], array.encoded->args[0])
          || (op == Store && !has_sort(args[2], array.encoded->args[1]))) {
        error(loc, refusal);
      }
      if (op == Store) {
        return typed(make_term(Store, terms, loc), array);
      }
      // A value outside the enumeration, which an array may hold, stands for
      // one of its values.
      return clamp(typed(make_term(Select, terms, loc), array.encoded->args[1]),
                   loc);
    }
    default: error(loc, refusal);
  }
}

Term Reader::make_equality(const Term & a,
                           const Term & b,
                           const MoxiSort & sort,
                           const string & name,
                           const Location & loc)
{
  const Enumeration * inexact = inexact_equality(sort);
  if (!inexact) {
    return make_term(Equal, { a, b }, loc);
  }
  // The solver's equality compares arrays indexed by an enumeration also at
  // the other values of the sort that encodes it, and tells apart elements
  // outside an enumeration that all stand for its first value. In MoXI,
  // arrays are equal if their elements at each index are, so arrays indexed
  // by an enumeration are compared at its values, one element at a time, as
  // select reads them.
  const MoxiSort & index = sort.encoded->args[0];
  const MoxiSort & element = sort.encoded->args[1];
  const Enumeration * indices =
      index.encoded ? index.encoded->enumeration : nullptr;
  if (!indices) {
    // Other indices are too many to compare one at a time.
    refuse_encoding(
        *inexact,
        name + " would tell apart equal arrays of sort " + sort.to_string(),
        loc);
  }
  TermVec equalities;
  for (const Term & value : indices->values) {
    const MoxiTerm a_element =
        clamp(typed(make_term(Select, { a, value }, loc), element), loc);
    const MoxiTerm b_element =
        clamp(typed(make_term(Select, { b, value }, loc), element), loc);
    equalities.push_back(
        make_equality(a_element.term, b_element.term, element, name, loc));
  }
  return equalities.size() == 1 ? equalities[0]
                                : make_term(And, equalities, loc);
}

MoxiTerm Reader::clamp(const MoxiTerm & term, const Location & loc) const
{
  const Term in_range = range(term.sort(), term.term);
  if (!in_range) {
    return term;
  }
  return { make_term(
               Ite,
               { in_range, term.term, term.encoded->enumeration->values[0] },
               loc),
           term.encoded };
}

bool Reader::is_formula(const MoxiTerm & term) const
{
  // A solver that aliases Bool and (_ BitVec 1) gives an enumeration of two
  // values that bit-vectors encode the sort Bool, too.
  return !term.encoded && term.term->get_sort() == bool_sort_;
}

bool Reader::has_sort(const MoxiTerm & term, const MoxiSort & sort) const
{
  if (term.encoded || sort.encoded) {
    return term.encoded == sort.encoded;
  }
  if (term.term->get_sort() == sort.sort) {
    return true;
  }
  // smt-switch's MathSAT and Yices2 backends give a real with an integer
  // value, e.g. 1.0, the sort Int, so an integer passes for a real with them.
  const SolverEnum se = solver_->get_solver_enum();
  return (se == MSAT || se == MSAT_INTERPOLATOR || se == YICES2)
         && sort.sort->get_sort_kind() == REAL
         && term.term->get_sort()->get_sort_kind() == INT;
}

Term Reader::make_number(const string & text,
                         const Sort & sort,
                         const Location & loc)
{
  try {
    const Term term = solver_->make_term(text, sort);
    numbers_.emplace(term, Number{ text, sort });
    return term;
  }
  catch (const SmtException & e) {
    error(loc, "cannot make number " + text + ": " + e.what());
  }
}

Term Reader::make_term(const Op & op,
                       const TermVec & args,
                       const Location & loc) const
{
  try {
    return solver_->make_term(op, args);
  }
  catch (const SmtException & e) {
    error(loc,
          "cannot apply " + op.to_string() + " to arguments of sort "
              + describe_sorts(args) + ": " + e.what());
  }
}

}  // namespace moxi
}  // namespace pono
