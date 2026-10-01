/*********************                                                        */
/*! \file moxi_encoder.cpp
** \verbatim
** Top contributors (to current version):
**   Po-Chun Chien
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** (in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Frontend for the Model Exchange Interlingua (MoXI).
**        See https://doi.org/10.1007/978-3-031-65627-9_10 and
**        https://github.com/ModelChecker/IL for more information.
**
**/

#include "frontends/moxi_encoder.h"

#include <cassert>
#include <unordered_map>

#include "frontends/moxi_reader.h"
#include "smt-switch/utils.h"
#include "utils/exceptions.h"
#include "utils/logger.h"
#include "utils/str_util.h"

using namespace smt;
using namespace std;

namespace pono {

/** @return the conjunction of the terms, true if there are none */
static Term conjunction(const SmtSolver & solver, const TermVec & terms)
{
  if (terms.empty()) {
    return solver->make_term(true);
  }
  return terms.size() == 1 ? terms[0] : solver->make_term(And, terms);
}

MoxiEncoder::MoxiEncoder(const string & filename,
                         RelationalTransitionSystem & rts,
                         size_t check_idx)
    : rts_(rts), reader_(make_unique<moxi::Reader>(filename, rts.solver()))
{
  const vector<moxi::Check> & checks = reader_->checks();
  if (checks.empty()) {
    throw PonoException("MoXI file " + filename
                        + " has no check-system command");
  }
  if (check_idx >= checks.size()) {
    throw PonoException("Check-system index " + to_string(check_idx)
                        + " is out of range: MoXI file " + filename + " has "
                        + to_string(checks.size()) + " check-system commands");
  }
  const moxi::Check & check = checks[check_idx];
  // MoXI asks the queries of a :queries attribute to be satisfiable with the
  // same values of the declared constants and functions, which checking them
  // one at a time cannot ensure.
  if (!check.queries_attributes.empty()) {
    reader_->error(check.queries_attributes.front(),
                   "the :queries attribute is not supported, as its queries"
                   " must share the values of the declared constants and"
                   " functions");
  }
  if (check.queries.empty()) {
    reader_->error(check.loc, "the check-system command has no query");
  }
  if (check.queries.size() > 1) {
    reader_->error(check.loc,
                   "the check-system command has "
                       + to_string(check.queries.size())
                       + " queries, but only one query per command is"
                         " supported");
  }
  query_name_ = check.queries[0].name;
  encode(check, check.queries[0]);
}

MoxiEncoder::~MoxiEncoder() = default;

void MoxiEncoder::encode(const moxi::Check & check, const moxi::Query & query)
{
  const SmtSolver & solver = rts_.solver();

  const moxi::CheckFormula * current = nullptr;
  vector<const moxi::CheckFormula *> assumptions;
  vector<const moxi::CheckFormula *> reachables;
  for (const size_t i : query.formulas) {
    const moxi::CheckFormula & formula = check.formulas[i];
    switch (formula.kind) {
      case moxi::FormulaKind::ASSUMPTION:
        assumptions.push_back(&formula);
        break;
      case moxi::FormulaKind::REACHABLE: reachables.push_back(&formula); break;
      case moxi::FormulaKind::CURRENT: current = &formula; break;
      case moxi::FormulaKind::FAIRNESS:
        reader_->error(formula.loc,
                       "query " + query.name + " has the fairness condition "
                           + formula.name
                           + ", but only queries without fairness conditions"
                             " are supported");
    }
  }

  // The variables of the check-system rename those of the system it checks,
  // which a declared constant, fixed over time, joins as a state that keeps
  // its value.
  vector<FlatVar *> actuals;
  for (const auto * group : { &check.inputs, &check.outputs, &check.locals }) {
    for (const moxi::Variable & var : *group) {
      vars_.push_back({ var.name, var.sort, var.curr, var.next, false });
      actuals.push_back(&vars_.back());
    }
  }
  vector<FlatVar *> constants;
  for (const moxi::Constant & constant : reader_->constants()) {
    vars_.push_back({ constant.name,
                      constant.sort,
                      constant.placeholder,
                      constant.placeholder,
                      true });
    constants.push_back(&vars_.back());
  }
  flatten(*check.system, actuals, "");

  // A variable has to be a state if a formula that holds in each state
  // refers to it, or if a formula refers to its next value. Otherwise only
  // the transitions refer to it, which lets it be an input.
  auto refer = [](FlatVar & var,
                  const Term & curr,
                  const Term & next,
                  const UnorderedTermSet & state_refs,
                  const UnorderedTermSet & step_refs) {
    var.in_state |= state_refs.count(curr) > 0;
    var.in_step |= step_refs.count(curr) > 0 || step_refs.count(next) > 0;
    var.next_in_step |= step_refs.count(next) > 0;
  };

  // The formulas of a system refer to its variables and to constants. Their
  // symbols are collected once, for all instances of the system, although
  // the solver can simplify some away once an instance's variables are
  // substituted, e.g. if it passes the same variable for two of the system's.
  // That variable then stays a state where it could have been an input.
  struct References
  {
    UnorderedTermSet state;
    UnorderedTermSet step;
  };
  unordered_map<const moxi::System *, References> system_refs;
  for (const Instance & instance : instances_) {
    const moxi::System & system = *instance.system;
    const auto [it, is_new] = system_refs.try_emplace(&system);
    References & refs = it->second;
    if (is_new) {
      // The query's initiality condition replaces the initial one.
      if (system.init && !current) {
        get_free_symbolic_consts(system.init, refs.state);
      }
      if (system.inv) {
        get_free_symbolic_consts(system.inv, refs.state);
      }
      if (system.trans) {
        get_free_symbolic_consts(system.trans, refs.step);
      }
    }
    for (size_t i = 0; i < instance.actuals.size(); ++i) {
      const moxi::Variable & formal = system.variable(i);
      refer(*instance.actuals[i],
            formal.curr,
            formal.next,
            refs.state,
            refs.step);
    }
    for (FlatVar * constant : constants) {
      refer(*constant,
            constant->curr_placeholder,
            constant->next_placeholder,
            refs.state,
            refs.step);
    }
  }

  // The formulas of the query refer to the variables of the check-system
  // and to constants.
  UnorderedTermSet next_values;
  for (const FlatVar * var : actuals) {
    next_values.insert(var->next_placeholder);
  }
  UnorderedTermSet state_refs;
  UnorderedTermSet step_refs;
  if (current) {
    get_free_symbolic_consts(current->term, state_refs);
  }
  vector<bool> is_step_assumption;
  for (const moxi::CheckFormula * assumption : assumptions) {
    UnorderedTermSet refs;
    get_free_symbolic_consts(assumption->term, refs);
    bool is_step = false;
    for (const Term & ref : refs) {
      is_step |= next_values.count(ref) > 0;
    }
    is_step_assumption.push_back(is_step);
    (is_step ? step_refs : state_refs).insert(refs.begin(), refs.end());
  }
  for (const moxi::CheckFormula * reachable : reachables) {
    UnorderedTermSet refs;
    get_free_symbolic_consts(reachable->term, refs);
    state_refs.insert(refs.begin(), refs.end());
    step_refs.insert(refs.begin(), refs.end());
  }
  for (const auto * group : { &actuals, &constants }) {
    for (FlatVar * var : *group) {
      refer(*var,
            var->curr_placeholder,
            var->next_placeholder,
            state_refs,
            step_refs);
    }
  }

  for (FlatVar & var : vars_) {
    if (var.frozen && !var.in_state && !var.in_step) {
      continue;
    }
    const bool is_state = var.frozen || var.in_state || var.next_in_step;
    var.curr = make_variable(var.name, var.sort, is_state);
    if (var.frozen) {
      rts_.assign_next(var.curr, var.curr);
    } else if (is_state) {
      var.next = rts_.next(var.curr);
    }
  }

  // Each instance contributes the formulas of its system, with the variables
  // of the transition system substituted for those of the system at once.
  TermVec init;
  TermVec trans;
  TermVec inv;
  for (const Instance & instance : instances_) {
    const moxi::System & system = *instance.system;
    UnorderedTermMap to_rts;
    for (size_t i = 0; i < instance.actuals.size(); ++i) {
      const moxi::Variable & formal = system.variable(i);
      const FlatVar & actual = *instance.actuals[i];
      to_rts[formal.curr] = actual.curr;
      if (actual.next) {
        to_rts[formal.next] = actual.next;
      }
    }
    for (const FlatVar * constant : constants) {
      if (constant->curr) {
        to_rts[constant->curr_placeholder] = constant->curr;
      }
    }
    if (system.init && !current) {
      init.push_back(solver->substitute(system.init, to_rts));
    }
    if (system.trans) {
      trans.push_back(solver->substitute(system.trans, to_rts));
    }
    if (system.inv) {
      inv.push_back(solver->substitute(system.inv, to_rts));
    }
  }

  UnorderedTermMap to_rts;
  for (const auto * group : { &actuals, &constants }) {
    for (const FlatVar * var : *group) {
      if (var->curr) {
        to_rts[var->curr_placeholder] = var->curr;
      }
      if (var->next) {
        to_rts[var->next_placeholder] = var->next;
      }
    }
  }
  auto encoded = [&](const Term & term) {
    return solver->substitute(term, to_rts);
  };
  if (current) {
    rts_.constrain_init(encoded(current->term));
  } else if (!init.empty()) {
    rts_.constrain_init(conjunction(solver, init));
  }
  if (!trans.empty()) {
    rts_.constrain_trans(conjunction(solver, trans));
  }
  if (!inv.empty()) {
    rts_.add_constraint(conjunction(solver, inv));
  }
  for (size_t i = 0; i < assumptions.size(); ++i) {
    const Term assumption = encoded(assumptions[i]->term);
    if (is_step_assumption[i]) {
      rts_.constrain_trans(assumption);
    } else {
      rts_.add_constraint(assumption);
    }
  }

  // With several reachability conditions, each may hold at a different step,
  // so a flag remembers for each one that it held before.
  TermVec reached;
  for (size_t i = 0; i < reachables.size(); ++i) {
    const Term condition = encoded(reachables[i]->term);
    if (reachables.size() == 1) {
      reached.push_back(condition);
      break;
    }
    const Term flag = rts_.make_generated_statevar("reached_" + to_string(i),
                                                   solver->make_sort(BOOL));
    rts_.constrain_init(solver->make_term(Not, flag));
    const Term held = solver->make_term(Or, flag, condition);
    rts_.constrain_trans(solver->make_term(Equal, rts_.next(flag), held));
    reached.push_back(held);
  }
  prop_ = solver->make_term(Not, conjunction(solver, reached));

  logger.log(1,
             "Encoded query {} of the check-system for system {}",
             query.name,
             check.system->name);
  logger.log(2,
             "The MoXI transition system has {} state and {} input variables",
             rts_.statevars().size(),
             rts_.inputvars().size());
}

void MoxiEncoder::flatten(const moxi::System & system,
                          const vector<FlatVar *> & actuals,
                          const string & prefix)
{
  assert(actuals.size() == system.num_variables());
  instances_.push_back({ &system, actuals });

  // Each instance of a subsystem has its own local variables.
  for (const moxi::Subsystem & subsystem : system.subsystems) {
    const string name = prefix + subsystem.name + ".";
    vector<FlatVar *> sub_actuals;
    for (const size_t pos : subsystem.args) {
      sub_actuals.push_back(actuals[pos]);
    }
    for (const moxi::Variable & local : subsystem.system->locals) {
      vars_.push_back(
          { name + local.name, local.sort, nullptr, nullptr, false });
      sub_actuals.push_back(&vars_.back());
    }
    flatten(*subsystem.system, sub_actuals, name);
  }
}

Term MoxiEncoder::make_variable(const string & name,
                                const Sort & sort,
                                bool is_state)
{
  // Names can clash, e.g. a check-system variable with a declared constant
  // it hides, or either with an uninterpreted function, and the solver keeps
  // symbols apart by their names.
  string candidate = name;
  for (size_t n = 1; taken_names_.count(candidate)
                     || reader_->function_names().count(candidate)
                     || rts_.named_terms().count(candidate);
       ++n) {
    candidate = name + "#" + to_string(n);
  }
  taken_names_.insert(candidate);
  return is_state ? rts_.make_statevar(candidate, sort)
                  : rts_.make_inputvar(candidate, sort);
}

}  // namespace pono
