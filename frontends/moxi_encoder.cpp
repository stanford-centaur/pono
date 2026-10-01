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
#include <cstdint>

#include "frontends/moxi_script.h"
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
                         size_t query_idx)
    : rts_(rts),
      script_(make_unique<moxi::Script>(filename, rts.solver())),
      query_idx_(query_idx)
{
  const moxi::Check * check = nullptr;
  const moxi::Query * query = nullptr;
  for (const moxi::Check & c : script_->checks()) {
    for (const moxi::Query & q : c.queries) {
      if (query_names_.size() == query_idx) {
        check = &c;
        query = &q;
      }
      query_names_.push_back(q.name);
    }
  }
  if (query_names_.empty()) {
    throw PonoException("MoXI file " + filename
                        + " has no check-system command with a query");
  }
  if (!query) {
    throw PonoException("Query index " + to_string(query_idx)
                        + " is out of range: MoXI file " + filename + " has "
                        + to_string(query_names_.size()) + " queries");
  }
  encode(*check, *query);
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
        script_->error(formula.loc,
                       "query " + query.name + " has the fairness condition "
                           + formula.name
                           + ", but only queries without fairness conditions"
                             " are supported");
    }
  }

  // The variables of the check-system rename those of the system it checks,
  // which a declared constant, fixed over time, joins as a state that keeps
  // its value.
  vector<const FlatVar *> actuals;
  for (const auto * group : { &check.inputs, &check.outputs, &check.locals }) {
    for (const moxi::Variable & var : *group) {
      vars_.push_back({ var.name, &var.sort, var.curr, var.next, false });
      actuals.push_back(&vars_.back());
    }
  }
  for (const moxi::Constant & constant : script_->constants()) {
    vars_.push_back({ constant.name,
                      &constant.sort,
                      constant.placeholder,
                      constant.placeholder,
                      true });
  }
  flatten(*check.system, actuals, "");

  // The query's initiality condition replaces the initial one.
  const Term init = current ? current->term : conjunction(solver, init_);
  const Term trans = conjunction(solver, trans_);
  const Term inv = conjunction(solver, inv_);

  // A variable has to be a state if a formula that holds in each state
  // refers to it, or if a formula refers to its next value. Otherwise only
  // the transitions refer to it, which lets it be an input.
  UnorderedTermSet next_values;
  for (const FlatVar & var : vars_) {
    if (!var.frozen) {
      next_values.insert(var.next);
    }
  }
  UnorderedTermSet state_refs;
  UnorderedTermSet step_refs;
  get_free_symbolic_consts(init, state_refs);
  get_free_symbolic_consts(inv, state_refs);
  get_free_symbolic_consts(trans, step_refs);
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

  UnorderedTermMap to_rts;
  for (const FlatVar & var : vars_) {
    const bool in_state = state_refs.count(var.curr) > 0;
    const bool in_step =
        step_refs.count(var.curr) > 0 || step_refs.count(var.next) > 0;
    if (var.frozen && !in_state && !in_step) {
      continue;
    }
    const bool is_state = var.frozen || in_state || step_refs.count(var.next);
    const Term term = make_variable(var.name, *var.sort, is_state);
    to_rts[var.curr] = term;
    if (var.frozen) {
      rts_.assign_next(term, term);
    } else if (is_state) {
      to_rts[var.next] = rts_.next(term);
    }
    if (const moxi::EnumSort * enumeration = var.sort->enumeration) {
      // The bit-vectors encoding an enumeration can outnumber its values.
      const uint64_t num_values = enumeration->values.size();
      const uint64_t width = var.sort->sort->get_width();
      if (width < 64 && num_values < (uint64_t{ 1 } << width)) {
        rts_.add_constraint(
            rts_.make_term(BVUle,
                           term,
                           rts_.make_term(static_cast<int64_t>(num_values - 1),
                                          var.sort->sort)));
      }
    }
  }

  auto encoded = [&](const Term & term) {
    return solver->substitute(term, to_rts);
  };
  if (current || !init_.empty()) {
    rts_.constrain_init(encoded(init));
  }
  if (!trans_.empty()) {
    rts_.constrain_trans(encoded(trans));
  }
  if (!inv_.empty()) {
    rts_.add_constraint(encoded(inv));
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
                          const vector<const FlatVar *> & actuals,
                          const string & prefix)
{
  assert(actuals.size() == system.num_variables());
  const SmtSolver & solver = rts_.solver();

  UnorderedTermMap substitution;
  for (size_t i = 0; i < actuals.size(); ++i) {
    const moxi::Variable & formal = system.variable(i);
    substitution[formal.curr] = actuals[i]->curr;
    substitution[formal.next] = actuals[i]->next;
  }
  if (system.init) {
    init_.push_back(solver->substitute(system.init, substitution));
  }
  if (system.trans) {
    trans_.push_back(solver->substitute(system.trans, substitution));
  }
  if (system.inv) {
    inv_.push_back(solver->substitute(system.inv, substitution));
  }

  // Each instance of a subsystem has its own local variables.
  for (const moxi::Subsystem & subsystem : system.subsystems) {
    const string name = prefix + subsystem.name + ".";
    vector<const FlatVar *> sub_actuals;
    for (const size_t pos : subsystem.args) {
      sub_actuals.push_back(actuals[pos]);
    }
    for (const moxi::Variable & local : subsystem.system->locals) {
      const Sort & sort = local.sort.sort;
      vars_.push_back({ name + local.name,
                        &local.sort,
                        script_->make_placeholder(sort),
                        script_->make_placeholder(sort),
                        false });
      sub_actuals.push_back(&vars_.back());
    }
    flatten(*subsystem.system, sub_actuals, name);
  }
}

Term MoxiEncoder::make_variable(const string & name,
                                const moxi::SortInfo & sort,
                                bool is_state)
{
  // Names can clash, e.g. a check-system variable with a declared constant
  // it hides, or either with an uninterpreted function, and the solver keeps
  // symbols apart by their names.
  string candidate = name;
  for (size_t n = 1; taken_names_.count(candidate)
                     || script_->function_names().count(candidate)
                     || rts_.named_terms().count(candidate);
       ++n) {
    candidate = name + "#" + to_string(n);
  }
  taken_names_.insert(candidate);
  return is_state ? rts_.make_statevar(candidate, sort.sort)
                  : rts_.make_inputvar(candidate, sort.sort);
}

}  // namespace pono
