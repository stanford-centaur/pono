/*********************                                                        */
/*! \file moxi_script.cpp
** \verbatim
** Top contributors (to current version):
**   Po-Chun Chien
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** (in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief The definitions a MoXI script makes, and the driver that the
**        Flex/Bison parser for MoXI fills them in through.
**
**/

#include "frontends/moxi_script.h"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cctype>
#include <limits>
#include <ostream>
#include <sstream>

#include "utils/exceptions.h"
#include "utils/str_util.h"

using namespace smt;
using namespace std;

namespace pono {
namespace moxi {

namespace {

/** How an operator that SMT-LIB defines as binary takes more arguments. */
enum class Fold
{
  NONE,      ///< it does not
  LEFT,      ///< :left-assoc, e.g. (- a b c) is (- (- a b) c)
  RIGHT,     ///< :right-assoc, e.g. (=> a b c) is (=> a (=> b c))
  CHAIN,     ///< :chainable, e.g. (< a b c) is (and (< a b) (< b c))
  PAIRWISE,  ///< :pairwise, e.g. (distinct a b c) holds if no two are equal
};

/** How an operator reconciles integer and real arguments. */
enum class Numeric
{
  NONE,   ///< it does not take numbers, or takes them as they are
  UNIFY,  ///< integers become reals if any argument is a real
  REAL,   ///< integers always become reals, as for real division
};

/** A predefined function that maps to an operator of the solver. */
struct Builtin
{
  PrimOp op;
  size_t min_args;
  size_t max_args;
  Fold fold;
  Numeric numeric;
};

constexpr size_t ANY = numeric_limits<size_t>::max();

const unordered_map<string, Builtin> & builtins()
{
  static const unordered_map<string, Builtin> table{
    // core
    { "not", { Not, 1, 1, Fold::NONE, Numeric::NONE } },
    { "and", { And, 1, ANY, Fold::LEFT, Numeric::NONE } },
    { "or", { Or, 1, ANY, Fold::LEFT, Numeric::NONE } },
    { "xor", { Xor, 2, ANY, Fold::LEFT, Numeric::NONE } },
    { "=>", { Implies, 2, ANY, Fold::RIGHT, Numeric::NONE } },
    { "=", { Equal, 2, ANY, Fold::CHAIN, Numeric::UNIFY } },
    { "distinct", { Distinct, 2, ANY, Fold::PAIRWISE, Numeric::UNIFY } },
    // MoXI adds a binary disequality
    { "!=", { Distinct, 2, 2, Fold::NONE, Numeric::UNIFY } },
    { "ite", { Ite, 3, 3, Fold::NONE, Numeric::UNIFY } },
    // integers and reals
    { "+", { Plus, 1, ANY, Fold::LEFT, Numeric::UNIFY } },
    { "-", { Minus, 1, ANY, Fold::LEFT, Numeric::UNIFY } },
    { "*", { Mult, 2, ANY, Fold::LEFT, Numeric::UNIFY } },
    { "/", { Div, 2, ANY, Fold::LEFT, Numeric::REAL } },
    { "div", { IntDiv, 2, ANY, Fold::LEFT, Numeric::NONE } },
    { "mod", { Mod, 2, 2, Fold::NONE, Numeric::NONE } },
    { "abs", { Abs, 1, 1, Fold::NONE, Numeric::NONE } },
    { "<", { Lt, 2, ANY, Fold::CHAIN, Numeric::UNIFY } },
    { "<=", { Le, 2, ANY, Fold::CHAIN, Numeric::UNIFY } },
    { ">", { Gt, 2, ANY, Fold::CHAIN, Numeric::UNIFY } },
    { ">=", { Ge, 2, ANY, Fold::CHAIN, Numeric::UNIFY } },
    { "to_real", { To_Real, 1, 1, Fold::NONE, Numeric::NONE } },
    { "to_int", { To_Int, 1, 1, Fold::NONE, Numeric::NONE } },
    { "is_int", { Is_Int, 1, 1, Fold::NONE, Numeric::NONE } },
    // bit-vectors
    { "concat", { Concat, 2, ANY, Fold::LEFT, Numeric::NONE } },
    { "bvnot", { BVNot, 1, 1, Fold::NONE, Numeric::NONE } },
    { "bvneg", { BVNeg, 1, 1, Fold::NONE, Numeric::NONE } },
    { "bvand", { BVAnd, 2, ANY, Fold::LEFT, Numeric::NONE } },
    { "bvor", { BVOr, 2, ANY, Fold::LEFT, Numeric::NONE } },
    { "bvxor", { BVXor, 2, ANY, Fold::LEFT, Numeric::NONE } },
    { "bvadd", { BVAdd, 2, ANY, Fold::LEFT, Numeric::NONE } },
    { "bvmul", { BVMul, 2, ANY, Fold::LEFT, Numeric::NONE } },
    { "bvnand", { BVNand, 2, 2, Fold::NONE, Numeric::NONE } },
    { "bvnor", { BVNor, 2, 2, Fold::NONE, Numeric::NONE } },
    { "bvxnor", { BVXnor, 2, 2, Fold::NONE, Numeric::NONE } },
    { "bvsub", { BVSub, 2, 2, Fold::NONE, Numeric::NONE } },
    { "bvudiv", { BVUdiv, 2, 2, Fold::NONE, Numeric::NONE } },
    { "bvsdiv", { BVSdiv, 2, 2, Fold::NONE, Numeric::NONE } },
    { "bvurem", { BVUrem, 2, 2, Fold::NONE, Numeric::NONE } },
    { "bvsrem", { BVSrem, 2, 2, Fold::NONE, Numeric::NONE } },
    { "bvsmod", { BVSmod, 2, 2, Fold::NONE, Numeric::NONE } },
    { "bvshl", { BVShl, 2, 2, Fold::NONE, Numeric::NONE } },
    { "bvlshr", { BVLshr, 2, 2, Fold::NONE, Numeric::NONE } },
    { "bvashr", { BVAshr, 2, 2, Fold::NONE, Numeric::NONE } },
    { "bvcomp", { BVComp, 2, 2, Fold::NONE, Numeric::NONE } },
    { "bvult", { BVUlt, 2, 2, Fold::NONE, Numeric::NONE } },
    { "bvule", { BVUle, 2, 2, Fold::NONE, Numeric::NONE } },
    { "bvugt", { BVUgt, 2, 2, Fold::NONE, Numeric::NONE } },
    { "bvuge", { BVUge, 2, 2, Fold::NONE, Numeric::NONE } },
    { "bvslt", { BVSlt, 2, 2, Fold::NONE, Numeric::NONE } },
    { "bvsle", { BVSle, 2, 2, Fold::NONE, Numeric::NONE } },
    { "bvsgt", { BVSgt, 2, 2, Fold::NONE, Numeric::NONE } },
    { "bvsge", { BVSge, 2, 2, Fold::NONE, Numeric::NONE } },
    { "ubv_to_int", { UBV_To_Int, 1, 1, Fold::NONE, Numeric::NONE } },
    { "sbv_to_int", { SBV_To_Int, 1, 1, Fold::NONE, Numeric::NONE } },
    { "bv2nat", { UBV_To_Int, 1, 1, Fold::NONE, Numeric::NONE } },
    // arrays
    { "select", { Select, 2, 2, Fold::NONE, Numeric::NONE } },
    { "store", { Store, 3, 3, Fold::NONE, Numeric::NONE } },
  };
  return table;
}

/** Predefined functions that the solver has no operator for. */
const unordered_set<string> derived_builtins{ "bvredand", "bvredor" };

/** Predefined functions that take indices, as in ((_ extract 7 0) x). */
const unordered_map<string, PrimOp> indexed_builtins{
  { "extract", Extract },
  { "zero_extend", Zero_Extend },
  { "sign_extend", Sign_Extend },
  { "rotate_left", Rotate_Left },
  { "rotate_right", Rotate_Right },
  { "repeat", Repeat },
  { "int2bv", Int_To_BV },
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
         || derived_builtins.count(name) || indexed_builtins.count(name);
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
bool parse_logic(const string & logic, bool & has_int, bool & has_real)
{
  if (logic == "ALL") {
    has_int = has_real = true;
    return true;
  }
  string rest = logic;
  if (syntax_analysis::StrStartsWith(rest, "QF_")) {
    rest = rest.substr(3);
  }
  bool has_theory = false;
  for (const string theory : { "AX", "A", "UF", "BV", "FP", "DT", "S" }) {
    if (syntax_analysis::StrStartsWith(rest, theory)) {
      rest = rest.substr(theory.size());
      has_theory = true;
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
  return true;
}

string describe(const SortInfo & sort)
{
  return sort.enumeration ? sort.enumeration->name : sort.sort->to_string();
}

/** @return whether a term of the actual sort can stand where the expected
 *  one is. An integer can stand for a real: some solvers, e.g. MathSAT, keep
 *  integral constants integers even when asked for reals, and mix the two. */
bool fits(const Sort & actual, const Sort & expected)
{
  return actual == expected
         || (expected->get_sort_kind() == REAL
             && actual->get_sort_kind() == INT);
}

/** @return a name for a placeholder symbol. A solver refuses a second symbol
 *  of a name, and scripts may share one, so the count is shared too. */
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

}  // namespace

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

Script::Script(const string & filename, const SmtSolver & solver)
    : filename_(filename), solver_(solver), bool_sort_(solver->make_sort(BOOL))
{
  parse();
}

void Script::error(const Location & loc, const string & message) const
{
  ostringstream text;
  text << filename_ << ":" << loc << ": " << message;
  throw PonoException(text.str());
}

// commands

void Script::set_logic(const string & logic, const Location & loc)
{
  if (logic_set_) {
    error(loc, "the logic is already set");
  }
  if (!parse_logic(logic, logic_has_int_, logic_has_real_)) {
    error(loc, "unknown logic " + logic);
  }
  logic_set_ = true;
}

void Script::declare_sort(const string & name,
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

void Script::define_sort(const string & name,
                         const vector<string> & params,
                         const SortExpr & body,
                         const Location & loc)
{
  check_new_sort_name(name, loc);
  // Resolve the body once to report its mistakes here, with the parameters
  // standing for an arbitrary sort.
  SortParams dummies;
  for (const string & param : params) {
    if (!dummies.emplace(param, SortInfo{ bool_sort_ }).second) {
      error(loc, "sort parameter " + param + " is repeated");
    }
  }
  resolve_sort(body, &dummies);
  defined_sorts_[name] = { params, body };
}

void Script::declare_enum_sort(const string & name,
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
  uint64_t width = 1;
  while (width < 64 && (uint64_t{ 1 } << width) < values.size()) {
    ++width;
  }
  enum_sorts_.push_back({ name, values, solver_->make_sort(BV, width) });
  const EnumSort & enumeration = enum_sorts_.back();
  enum_sort_names_[name] = &enumeration;
  for (size_t i = 0; i < values.size(); ++i) {
    enum_values_[values[i]] =
        solver_->make_term(static_cast<int64_t>(i), enumeration.sort);
  }
}

void Script::declare_fun(const string & name,
                         const vector<SortExpr> & args,
                         const SortExpr & result,
                         const Location & loc)
{
  check_new_global_name(name, loc);
  const SortInfo result_sort = resolve_sort(result);
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

void Script::begin_define_fun(const string & name,
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
    const SortInfo sort = resolve_sort(param.sort);
    macro_.params.push_back(make_placeholder(sort.sort));
    macro_.param_sorts.push_back(sort);
    bind(param.name, macro_.params.back());
  }
  macro_.result = resolve_sort(result);
  next_allowed_ = false;
  formula_attribute_ = "define-fun";
}

void Script::end_define_fun(const Term & body, const Location & loc)
{
  pop_scope();
  macro_.body = coerce(body, macro_.result.sort, loc);
  if (!fits(macro_.body->get_sort(), macro_.result.sort)) {
    error(loc,
          "the definition of " + macro_name_ + " has sort "
              + body->get_sort()->to_string() + ", but "
              + describe(macro_.result) + " is declared");
  }
  macros_[macro_name_] = std::move(macro_);
}

void Script::begin_system(const string & name, const Location & loc)
{
  scope_ = Scope::SYSTEM;
  system_ = make_shared<System>();
  system_->name = name;
  system_->loc = loc;
  declared_attributes_.clear();
  command_names_.clear();
  variables_in_scope_.clear();
}

void Script::end_system()
{
  // A later definition of the same name replaces this one, while the
  // systems defined in between keep instantiating this one.
  systems_[system_->name] = system_;
  system_.reset();
  variables_in_scope_.clear();
  scope_ = Scope::TOP;
}

void Script::begin_check(const string & system, const Location & loc)
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

void Script::end_check()
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

vector<Variable> & Script::variables(const string & attribute,
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

void Script::declare_variables(const string & attribute,
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
    const SortInfo sort = resolve_sort(var.sort);
    group.push_back({ var.name,
                      sort,
                      make_placeholder(sort.sort),
                      make_placeholder(sort.sort),
                      var.loc });
  }
}

void Script::check_signature(const string & attribute,
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
                + describe(declared[i].sort) + ", but the variable "
                + expected[i].name + " of system " + system.name
                + " it renames has sort " + describe(expected[i].sort));
    }
  }
}

void Script::enter_variables(const vector<Variable> & vars)
{
  for (const Variable & var : vars) {
    variables_in_scope_[var.name] = &var;
  }
}

void Script::end_variables()
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

void Script::begin_formula(const string & attribute)
{
  next_allowed_ = attribute == ":trans" || attribute == ":assumption"
                  || attribute == ":fairness" || attribute == ":reachable";
  formula_attribute_ = attribute;
}

void Script::set_system_formula(const string & attribute,
                                const Term & formula,
                                const Location & loc)
{
  next_allowed_ = false;
  if (!declared_attributes_.insert(attribute).second) {
    error(loc, "attribute " + attribute + " appears more than once");
  }
  if (formula->get_sort()->get_sort_kind() != BOOL) {
    error(loc,
          "the formula of " + attribute + " has sort "
              + formula->get_sort()->to_string() + " instead of Bool");
  }
  if (attribute == ":init") {
    system_->init = formula;
  } else if (attribute == ":trans") {
    system_->trans = formula;
  } else {
    assert(attribute == ":inv");
    system_->inv = formula;
  }
}

void Script::add_subsystem(const string & name,
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
                + describe(actual.sort) + " for " + formal.name + " of sort "
                + describe(formal.sort));
    }
    subsystem.args.push_back(pos->second);
  }
  system_->subsystems.push_back(std::move(subsystem));
}

void Script::add_check_formula(const string & attribute,
                               const string & name,
                               const Term & formula,
                               const Location & loc)
{
  next_allowed_ = false;
  if (formula->get_sort()->get_sort_kind() != BOOL) {
    error(loc,
          "the formula " + name + " has sort "
              + formula->get_sort()->to_string() + " instead of Bool");
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
  check_.formulas.push_back({ name, kind, formula, loc });
}

void Script::add_query(const QueryExpr & query)
{
  // Queries may name formulas that come later, so they are resolved at the
  // end of the command.
  query_exprs_.push_back(query);
}

// sorts

Sort Script::int_sort(const Location & loc)
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

Sort Script::real_sort(const Location & loc)
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

Sort Script::numeral_sort(const Location & loc)
{
  // In SMT-LIB, numerals are reals in the logics without integers.
  return logic_has_real_ && !logic_has_int_ ? real_sort(loc) : int_sort(loc);
}

uint64_t Script::parse_index(const string & index, const Location & loc) const
{
  if (index.empty() || index.size() > 19
      || !all_of(index.begin(), index.end(), [](char c) {
           return isdigit(static_cast<unsigned char>(c));
         })) {
    error(loc, "expected a numeral that fits 64 bits, not " + index);
  }
  return stoull(index);
}

SortInfo Script::resolve_sort(const SortExpr & expr, const SortParams * params)
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
      return { bool_sort_ };
    } else if (name == "Int" && plain) {
      return { int_sort(expr.loc) };
    } else if (name == "Real" && plain) {
      return { real_sort(expr.loc) };
    } else if (name == "BitVec" && expr.indices.size() == 1
               && expr.args.empty()) {
      const uint64_t width = parse_index(expr.indices[0], expr.loc);
      if (!width) {
        error(expr.loc, "bit-vectors need a positive width");
      }
      return { solver_->make_sort(BV, width) };
    } else if (name == "Array" && expr.indices.empty()
               && expr.args.size() == 2) {
      return { solver_->make_sort(ARRAY,
                                  resolve_sort(expr.args[0], params).sort,
                                  resolve_sort(expr.args[1], params).sort) };
    }

    auto enumeration = enum_sort_names_.find(name);
    if (enumeration != enum_sort_names_.end() && plain) {
      return { enumeration->second->sort, enumeration->second };
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
        return { declared->second.sort };
      }
      SortVec args;
      for (const SortExpr & arg : expr.args) {
        args.push_back(resolve_sort(arg, params).sort);
      }
      return { solver_->make_sort(declared->second.sort, args) };
    }
  }
  catch (const SmtException & e) {
    error(expr.loc, "cannot make sort " + name + ": " + e.what());
  }

  const bool known = builtin_sorts.count(name) || enum_sort_names_.count(name)
                     || defined_sorts_.count(name)
                     || declared_sorts_.count(name);
  error(expr.loc,
        known ? "wrong indices or arguments for sort " + name
              : "unknown sort " + name);
}

void Script::check_new_sort_name(const string & name,
                                 const Location & loc) const
{
  if (builtin_sorts.count(name) || enum_sort_names_.count(name)
      || defined_sorts_.count(name) || declared_sorts_.count(name)) {
    error(loc, "sort " + name + " is already defined");
  }
}

void Script::check_new_global_name(const string & name,
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

Term Script::make_placeholder(const Sort & sort)
{
  return solver_->make_symbol(placeholder_name(), sort);
}

// terms

const Variable * Script::lookup_variable(const string & name) const
{
  auto it = variables_in_scope_.find(name);
  return it == variables_in_scope_.end() ? nullptr : it->second;
}

Term Script::lookup_bound(const string & name) const
{
  auto it = bound_.find(name);
  return it == bound_.end() ? Term() : it->second.back();
}

void Script::bind(const string & name, const Term & term)
{
  bound_[name].push_back(term);
  bound_scopes_.back().push_back(name);
}

void Script::push_let(const vector<VarBinding> & bindings)
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

void Script::pop_scope()
{
  for (const string & name : bound_scopes_.back()) {
    // bind() put the name in, so this does not insert it
    vector<Term> & terms = bound_[name];
    assert(!terms.empty());
    if (terms.size() == 1) {
      bound_.erase(name);
    } else {
      terms.pop_back();
    }
  }
  bound_scopes_.pop_back();
}

void Script::push_quantifier(const vector<SortedVar> & vars)
{
  bound_scopes_.emplace_back();
  TermVec params;
  unordered_set<string> seen;
  for (const SortedVar & var : vars) {
    if (!seen.insert(var.name).second) {
      error(var.loc, "quantifier binds " + var.name + " more than once");
    }
    const Sort sort = resolve_sort(var.sort).sort;
    params.push_back(solver_->make_param(placeholder_name(), sort));
    bind(var.name, params.back());
  }
  quantifier_params_.push_back(std::move(params));
}

Term Script::pop_quantifier(bool forall,
                            const Term & body,
                            const Location & loc)
{
  TermVec args = std::move(quantifier_params_.back());
  quantifier_params_.pop_back();
  pop_scope();
  args.push_back(body);
  return make_term(forall ? Forall : Exists, args, loc);
}

Term Script::make_numeral(const string & digits, const Location & loc)
{
  try {
    return solver_->make_term(digits, numeral_sort(loc));
  }
  catch (const SmtException & e) {
    error(loc, "cannot make numeral " + digits + ": " + e.what());
  }
}

Term Script::make_decimal(const string & text, const Location & loc)
{
  // Solvers take a decimal as a fraction more readily than as is.
  const size_t point = text.find('.');
  string denominator = "1" + string(text.size() - point - 1, '0');
  string numerator = text.substr(0, point) + text.substr(point + 1);
  numerator.erase(0,
                  min(numerator.find_first_not_of('0'), numerator.size() - 1));
  try {
    return solver_->make_term(numerator + "/" + denominator, real_sort(loc));
  }
  catch (const SmtException & e) {
    error(loc, "cannot make decimal " + text + ": " + e.what());
  }
}

Term Script::make_binary(const string & bits, const Location & loc)
{
  try {
    return solver_->make_term(bits, solver_->make_sort(BV, bits.size()), 2);
  }
  catch (const SmtException & e) {
    error(loc, "cannot make bit-vector #b" + bits + ": " + e.what());
  }
}

Term Script::make_hexadecimal(const string & digits, const Location & loc)
{
  try {
    return solver_->make_term(
        digits, solver_->make_sort(BV, 4 * digits.size()), 16);
  }
  catch (const SmtException & e) {
    error(loc, "cannot make bit-vector #x" + digits + ": " + e.what());
  }
}

Term Script::make_next(const string & name, const Location & loc)
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
  return var->next;
}

Term Script::make_identifier_term(const Identifier & id)
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
    const SortInfo sort = resolve_sort(*id.qualifier);
    const Term term = coerce(make_identifier_term(plain), sort.sort, id.loc);
    if (!fits(term->get_sort(), sort.sort)) {
      error(id.loc, id.name + " does not have sort " + describe(sort));
    }
    return term;
  }
  if (id.primed) {
    return make_next(id.name, id.loc);
  }

  const string & name = id.name;
  if (name == "true" || name == "false") {
    return solver_->make_term(name == "true");
  }
  if (Term term = lookup_bound(name)) {
    return term;
  }
  if (const Variable * var = lookup_variable(name)) {
    return var->curr;
  }
  auto value = enum_values_.find(name);
  if (value != enum_values_.end()) {
    return value->second;
  }
  auto constant = constant_names_.find(name);
  if (constant != constant_names_.end()) {
    return constant->second->placeholder;
  }
  auto macro = macros_.find(name);
  if (macro != macros_.end()) {
    if (!macro->second.params.empty()) {
      error(id.loc, "function " + name + " needs arguments");
    }
    return macro->second.body;
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

Term Script::make_indexed_constant(const Identifier & id)
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
    return solver_->make_term(id.name.substr(2), solver_->make_sort(BV, width));
  }
  catch (const SmtException & e) {
    error(id.loc, "cannot make bit-vector " + id.name + ": " + e.what());
  }
}

Term Script::make_application(const Identifier & id,
                              const TermVec & args,
                              const Location & loc)
{
  if (id.primed) {
    error(id.loc, "the primed symbol " + id.name + "' is not a function");
  }
  if (id.qualifier) {
    const SortInfo sort = resolve_sort(*id.qualifier);
    if (id.name == "const" && id.indices.empty()) {
      if (sort.sort->get_sort_kind() != ARRAY) {
        error(id.loc, "(as const S) needs an array sort S");
      }
      if (args.size() != 1) {
        error(loc, "a constant array takes a single value");
      }
      const Term elem = coerce(args[0], sort.sort->get_elemsort(), loc);
      try {
        return solver_->make_term(elem, sort.sort);
      }
      catch (const SmtException & e) {
        error(loc, string("cannot make constant array: ") + e.what());
      }
    }
    // Otherwise the sort only disambiguates, which the arguments do already.
    Identifier plain = id;
    plain.qualifier.reset();
    const Term term =
        coerce(make_application(plain, args, loc), sort.sort, loc);
    if (!fits(term->get_sort(), sort.sort)) {
      error(loc,
            "the application of " + id.name + " does not have sort "
                + describe(sort));
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
      children.push_back(coerce(args[i], f.arg_sorts[i].sort, loc));
    }
    return make_term(Apply, children, loc);
  }
  return make_builtin(id.name, args, loc);
}

Term Script::apply_macro(const string & name,
                         const Macro & macro,
                         TermVec args,
                         const Location & loc)
{
  if (args.size() != macro.params.size()) {
    error(loc,
          "function " + name + " takes " + to_string(macro.params.size())
              + " arguments, not " + to_string(args.size()));
  }
  UnorderedTermMap substitution;
  for (size_t i = 0; i < args.size(); ++i) {
    const Term arg = coerce(args[i], macro.param_sorts[i].sort, loc);
    if (!fits(arg->get_sort(), macro.param_sorts[i].sort)) {
      error(loc,
            "argument " + to_string(i + 1) + " of " + name + " has sort "
                + arg->get_sort()->to_string() + " instead of "
                + describe(macro.param_sorts[i]));
    }
    substitution[macro.params[i]] = arg;
  }
  try {
    return solver_->substitute(macro.body, substitution);
  }
  catch (const SmtException & e) {
    error(loc, "cannot apply " + name + ": " + e.what());
  }
}

Term Script::make_indexed_application(const Identifier & id,
                                      const TermVec & args,
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
  if (id.name == "divisible") {
    if (!indices[0]) {
      error(id.loc, "divisible needs a positive index");
    }
    const Term divisor =
        solver_->make_term(to_string(indices[0]), int_sort(id.loc));
    const Term zero = solver_->make_term(int64_t{ 0 }, int_sort(id.loc));
    return make_term(
        Equal, { make_term(Mod, { args[0], divisor }, loc), zero }, loc);
  }
  const Op op = num_indices == 2 ? Op(it->second, indices[0], indices[1])
                                 : Op(it->second, indices[0]);
  return make_term(op, args, loc);
}

Term Script::make_builtin(const string & name,
                          TermVec args,
                          const Location & loc)
{
  if (derived_builtins.count(name)) {
    if (args.size() != 1 || args[0]->get_sort()->get_sort_kind() != BV) {
      error(loc, name + " takes a single bit-vector");
    }
    const Sort sort = args[0]->get_sort();
    const Sort bit = solver_->make_sort(BV, 1);
    // bvredand: all bits are set; bvredor: some bit is set
    const Term zero = solver_->make_term(int64_t{ 0 }, sort);
    const Term holds =
        name == "bvredand"
            ? make_term(
                  Equal, { args[0], make_term(BVNot, { zero }, loc) }, loc)
            : make_term(Distinct, { args[0], zero }, loc);
    return make_term(Ite,
                     { holds,
                       solver_->make_term(int64_t{ 1 }, bit),
                       solver_->make_term(int64_t{ 0 }, bit) },
                     loc);
  }

  auto it = builtins().find(name);
  if (it == builtins().end()) {
    error(loc, "unknown function " + name);
  }
  const Builtin & builtin = it->second;
  if (args.size() < builtin.min_args || args.size() > builtin.max_args) {
    error(loc,
          name + " cannot take " + to_string(args.size()) + " argument"
              + (args.size() == 1 ? "" : "s"));
  }

  if (builtin.numeric == Numeric::REAL) {
    for (Term & arg : args) {
      arg = coerce(arg, real_sort(loc), loc);
    }
  } else if (builtin.numeric == Numeric::UNIFY) {
    unify_numeric(args, loc);
  }

  if (builtin.op == Minus && args.size() == 1) {
    return make_term(Negate, args, loc);
  }
  if (args.size() == 1 && builtin.fold == Fold::LEFT) {
    return args[0];
  }
  switch (builtin.fold) {
    case Fold::NONE: return make_term(builtin.op, args, loc);
    case Fold::LEFT: {
      if (is_variadic(builtin.op)) {
        return make_term(builtin.op, args, loc);
      }
      Term result = args[0];
      for (size_t i = 1; i < args.size(); ++i) {
        result = make_term(builtin.op, { result, args[i] }, loc);
      }
      return result;
    }
    case Fold::RIGHT: {
      Term result = args.back();
      for (size_t i = args.size() - 1; i-- > 0;) {
        result = make_term(builtin.op, { args[i], result }, loc);
      }
      return result;
    }
    case Fold::CHAIN:
    case Fold::PAIRWISE: {
      TermVec conjuncts;
      for (size_t i = 0; i + 1 < args.size(); ++i) {
        const size_t last =
            builtin.fold == Fold::CHAIN ? i + 1 : args.size() - 1;
        for (size_t j = i + 1; j <= last; ++j) {
          conjuncts.push_back(make_term(builtin.op, { args[i], args[j] }, loc));
        }
      }
      return conjuncts.size() == 1 ? conjuncts[0]
                                   : make_term(And, conjuncts, loc);
    }
  }
  assert(false);
  return Term();
}

void Script::unify_numeric(TermVec & args, const Location & loc) const
{
  bool has_real = false;
  for (const Term & arg : args) {
    has_real |= arg->get_sort()->get_sort_kind() == REAL;
  }
  if (!has_real) {
    return;
  }
  for (Term & arg : args) {
    if (arg->get_sort()->get_sort_kind() == INT) {
      arg = make_term(To_Real, { arg }, loc);
    }
  }
}

Term Script::coerce(const Term & term,
                    const Sort & sort,
                    const Location & loc) const
{
  if (sort->get_sort_kind() == REAL
      && term->get_sort()->get_sort_kind() == INT) {
    return make_term(To_Real, { term }, loc);
  }
  return term;
}

Term Script::make_term(const Op & op,
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
