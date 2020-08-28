/*********************                                                        */
/*! \file ceg_prophecy.cpp
** \verbatim
** Top contributors (to current version):
**   Ahmed Irfan, Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief An implementation of Counter-Example Guided Prophecy for array
**        model checking. It is parameterized by an underlying model checking
**        procedure which need not handle arrays (only UF). However, a common
**        instantiation is with an IC3-style procedure, in which case we
**        often refer to this algorithm as "prophic3".
**
**/

#include "assert.h"

#include "engines/ceg_prophecy.h"
#include "utils/make_provers.h"

using namespace smt;
using namespace std;

namespace pono {

CegProphecy::CegProphecy(const Property & p, Engine e, smt::SolverEnum se)
    : super(p, se),
      conc_ts_(p.transition_system()),
      solver_(conc_ts_.solver()),
      abs_ts_(solver_),
      e_(e),
      abs_unroller_(abs_ts_, solver_),
      aa_(conc_ts_, abs_ts_, true),
      aae_(p, aa_, abs_unroller_),
      pm_(abs_ts_)
{
  initialize();
}

CegProphecy::CegProphecy(const Property & p, Engine e, const SmtSolver & solver)
    : super(p, solver),
      conc_ts_(p.transition_system()),
      solver_(conc_ts_.solver()),
      abs_ts_(solver_),
      e_(e),
      abs_unroller_(abs_ts_, solver_),
      aa_(conc_ts_, abs_ts_, true),
      aae_(p, aa_, abs_unroller_),
      pm_(abs_ts_)
{
  initialize();
}

CegProphecy::CegProphecy(const PonoOptions & opt,
                         const Property & p,
                         Engine e,
                         smt::SolverEnum se)
    : super(opt, p, se),
      conc_ts_(p.transition_system()),
      solver_(conc_ts_.solver()),
      abs_ts_(solver_),
      e_(e),
      abs_unroller_(abs_ts_, solver_),
      aa_(conc_ts_, abs_ts_, true),
      aae_(p, aa_, abs_unroller_),
      pm_(abs_ts_)
{
  initialize();
}

CegProphecy::CegProphecy(const PonoOptions & opt,
                         const Property & p,
                         Engine e,
                         const smt::SmtSolver & solver)
    : super(opt, p, solver),
      conc_ts_(p.transition_system()),
      solver_(conc_ts_.solver()),
      abs_ts_(solver_),
      abs_unroller_(abs_ts_, solver_),
      e_(e),
      aa_(conc_ts_, abs_ts_, true),
      aae_(p, aa_, abs_unroller_),
      pm_(abs_ts_)
{
  initialize();
}

ProverResult CegProphecy::prove()
{
  ProverResult res = ProverResult::FALSE;
  while (res == ProverResult::FALSE) {
    // Refine the system
    // heuristic -- stop refining when no new axioms are needed.
    do {
      if (!refine()) {
        return ProverResult::FALSE;
      }
      reached_k_++;
    } while (num_added_axioms_);

    Property latest_prop(abs_ts_, solver_->make_term(Not, bad_));
    shared_ptr<Prover> prover =
        make_prover(e_, latest_prop, solver_->get_solver_enum(), options_);
    res = prover->prove();
  }

  return res;
}

ProverResult CegProphecy::check_until(int k)
{
  ProverResult res = ProverResult::FALSE;
  while (res == ProverResult::FALSE && reached_k_ <= k) {
    // Refine the system
    // heuristic -- stop refining when no new axioms are needed.
    do {
      if (!refine()) {
        return ProverResult::FALSE;
      }
      reached_k_++;
    } while (num_added_axioms_ && reached_k_ <= k);

    Property latest_prop(abs_ts_, solver_->make_term(Not, bad_));
    shared_ptr<Prover> prover =
        make_prover(e_, latest_prop, solver_->get_solver_enum(), options_);

    res = prover->check_until(k);
  }
  return res;
}

void CegProphecy::initialize()
{
  super::initialize();
  abstract();
}

void CegProphecy::abstract()
{
  // the ArrayAbstractor already abstracted the transition system on
  // construction -- only need to abstract bad
  bad_ = aa_.abstract(bad_);

  // the abstract system should have all the same state and input variables
  // but abstracted
  // plus it will have some new variables for the witnesses and lambdas
  assert(abs_ts_.statevars().size() >= conc_ts_.statevars().size());
  assert(abs_ts_.inputvars().size() >= conc_ts_.inputvars().size());
}

bool CegProphecy::refine()
{
  num_added_axioms_ = 0;
  // TODO use ArrayAxiomEnumerator and modifiers to refine the system
  // create BMC formula
  Term abs_bmc_formula = abs_unroller_.at_time(abs_ts_.init(), 0);
  for (int k = 0; k <= reached_k_; ++k) {
    abs_bmc_formula = solver_->make_term(
        And, abs_bmc_formula, abs_unroller_.at_time(abs_ts_.trans(), k));
  }
  abs_bmc_formula = solver_->make_term(
      And, abs_bmc_formula, abs_unroller_.at_time(bad_, reached_k_ + 1));

  // check array axioms over the abstract system
  if (!aae_.enumerate_axioms(abs_bmc_formula, reached_k_ + 1)) {
    // concrete CEX
    return false;
  }

  TermVec consecutive_axioms = aae_.get_consecutive_axioms();
  const AxiomVec & nonconsecutive_axioms = aae_.get_nonconsecutive_axioms();

  bool found_nonconsecutive_axioms = nonconsecutive_axioms.size();

  // TODO reduce axioms with an unsat core or specialized dropping

  if (found_nonconsecutive_axioms) {
    // TODO
    // prophecize and update ts
    // then look for axioms again

    // TODO: consider adding these axioms directly over the prophecy
    //       variables at the correct time
    //       for now, easier to just search for consecutive axioms

    // vector of pairs
    // first: prophecy variable
    // second: target (a history variable for non-zero delay)
    vector<pair<Term, Term>> proph_vars;
    for (AxiomInstantiation ax_inst : nonconsecutive_axioms) {
      assert(ax_inst.instantiations.size()
             == 1);  // expecting only one index instantiated
      for (auto idx : ax_inst.instantiations) {
        // number of steps before the property violation
        size_t delay = reached_k_ + 1 - abs_unroller_.get_curr_time(idx);
        // Prophecy Modifier will add prophecy and history variables
        // automatically here but it does NOT update the property
        proph_vars.push_back(pm_.get_proph(idx, delay));
      }
    }

    // now update bad_
    // the property would be updated as
    // (proph1=target1 /\ ... /\ prophn=targetn) -> prop
    // but because we're working with bad, this is equivalent to
    // (proph1=target1 /\ ... /\ prophn=targetn) /\ bad
    // where bad = !prop
    for (auto p : proph_vars) {
      Term proph_var = p.first;
      Term target = p.second;
      bad_ = solver_->make_term(
          And, solver_->make_term(Equal, proph_var, target), bad_);
    }

    // search for axioms again but don't include nonconsecutive ones
    bool ok = aae_.enumerate_axioms(abs_bmc_formula, reached_k_ + 1, false);
    // should be guaranteed to rule out counterexamples at this bound
    assert(ok);
    consecutive_axioms = aae_.get_consecutive_axioms();
    assert(!aae_.get_nonconsecutive_axioms().size());
  }

  // add consecutive axioms to the system
  // TODO: make sure we're adding current / next correctly
  for (auto ax : consecutive_axioms) {
    num_added_axioms_++;
    if (reached_k_ == -1) {
      // if only checking initial state
      // need to add to init
      abs_ts_.constrain_init(ax);
    }

    abs_ts_.constrain_trans(ax);
    if (abs_ts_.only_curr(ax)) {
      // add the next state version if it's an invariant over current state vars
      abs_ts_.constrain_trans(abs_ts_.next(ax));
    }
  }

  // TODO: update abs_ts_
  // TODO: set num_added_axioms_
  // TODO: update bad_

  // able to successfully refine
  return true;
}

}  // namespace pono
