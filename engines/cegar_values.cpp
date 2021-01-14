/*********************                                                        */
/*! \file cegar_values.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief A simple CEGAR loop that abstracts values with frozen variables
**        and refines by constraining the variable to the value again
**
**/

#include "engines/cegar_values.h"

#include <unordered_set>

#include "core/fts.h"
#include "core/rts.h"
#include "engines/ceg_prophecy_arrays.h"
#include "engines/ic3ia.h"
#include "smt-switch/identity_walker.h"
#include "smt/available_solvers.h"
#include "utils/exceptions.h"
#include "utils/logger.h"
#include "utils/make_provers.h"

using namespace smt;
using namespace std;

namespace pono {

// TODO add a value abstractor
//      make sure not to introduce nonlinearities
//      implement generic backend
//      then specialize for IC3IA

static unordered_set<PrimOp> nl_ops({ Mult,
                                      Div,
                                      Mod,
                                      Abs,
                                      Pow,
                                      IntDiv,
                                      BVMul,
                                      BVUdiv,
                                      BVSdiv,
                                      BVUrem,
                                      BVSrem,
                                      BVSmod });

/** Class for abstracting values with frozen variables
 */
class ValueAbstractor : public smt::IdentityWalker
{
 public:
  ValueAbstractor(TransitionSystem & ts, UnorderedTermMap & abstracted_values)
      : smt::IdentityWalker(ts.solver(), false),
        ts_(ts),
        abstracted_values_(abstracted_values)
  {
  }

 protected:
  smt::WalkerStepResult visit_term(smt::Term & term) override
  {
    if (!preorder_) {
      if (term->is_value() && term->get_sort()->get_sort_kind() != ARRAY) {
        // create a frozen variable
        Term frozen_var =
            ts_.make_statevar("__abs_" + term->to_string(), term->get_sort());
        // will be made frozen later
        // since the transition system will be modified after this
        save_in_cache(term, frozen_var);
        abstracted_values_[frozen_var] = term;
        return Walker_Continue;
      }

      Op op = term->get_op();
      if (!op.is_null() && nl_ops.find(op.prim_op) == nl_ops.end()) {
        // only rebuild terms with operators that can't
        // create nonlinearities by replacing constant values
        TermVec cached_children;
        Term c;
        for (const auto & t : term) {
          c = t;
          query_cache(t, c);
          cached_children.push_back(c);
        }
        save_in_cache(term, solver_->make_term(op, cached_children));
      } else {
        save_in_cache(term, term);
      }
    }
    return Walker_Continue;
  }

  TransitionSystem & ts_;
  UnorderedTermMap & abstracted_values_;
};

TransitionSystem create_fresh_ts(bool functional, const SmtSolver & solver)
{
  if (functional) {
    return FunctionalTransitionSystem(solver);
  } else {
    return RelationalTransitionSystem(solver);
  }
}

template <class Prover_T>
CegarValues<Prover_T>::CegarValues(const Property & p,
                                   const TransitionSystem & ts,
                                   const smt::SmtSolver & solver,
                                   PonoOptions opt)
    : super(p, create_fresh_ts(ts.is_functional(), solver), solver, opt),
      conc_ts_(ts, super::to_prover_solver_),
      prover_ts_(super::prover_interface_ts()),
      cegval_solver_(create_solver(solver->get_solver_enum())),
      to_cegval_solver_(cegval_solver_),
      from_cegval_solver_(super::solver_),
      cegval_ts_(cegval_solver_),
      cegval_un_(cegval_ts_)
{
  cegval_solver_->set_opt("produce-unsat-cores", "true");
}

template <class Prover_T>
ProverResult CegarValues<Prover_T>::check_until(int k)
{
  initialize();

  ProverResult res = ProverResult::FALSE;
  while (res == ProverResult::FALSE) {
    // need to call parent's check_until in case it
    // is another cegar loop rather than an engine
    res = super::check_until(k);

    if (res == ProverResult::FALSE) {
      if (!cegar_refine()) {
        return ProverResult::FALSE;
      }
    }
  }

  if (res == ProverResult::TRUE) {
    // update the invariant
    UnorderedTermMap super_to_vals;
    for (const auto & elem : to_vals_) {
      super_to_vals[from_cegval_solver_.transfer_term(elem.first)] =
          from_cegval_solver_.transfer_term(elem.second);
    }
    super::invar_ = super::solver_->substitute(super::invar_, super_to_vals);
  }

  return res;
}

template <class Prover_T>
void CegarValues<Prover_T>::initialize()
{
  if (super::initialized_) {
    return;
  }

  // specify which cegar_abstract in case
  // we're inheriting from another cegar algorithm
  CegarValues::cegar_abstract();
  super::initialize();

  // update local version of ts over fresh solver
  cegval_ts_ = TransitionSystem(prover_ts_, to_cegval_solver_);
  // update cache
  UnorderedTermMap & cache = from_cegval_solver_.get_cache();
  Term nv;
  for (const auto & sv : prover_ts_.statevars()) {
    nv = prover_ts_.next(sv);
    cache[to_cegval_solver_.transfer_term(sv)] = sv;
    cache[to_cegval_solver_.transfer_term(nv)] = nv;
  }

  for (const auto & iv : prover_ts_.inputvars()) {
    cache[to_cegval_solver_.transfer_term(iv)] = iv;
  }

  // create labels for each abstract value
  Sort boolsort = cegval_solver_->make_sort(BOOL);
  Term lbl;
  for (const auto & elem : to_vals_) {
    lbl = cegval_solver_->make_symbol("__assump_" + elem.second->to_string(),
                                      boolsort);
    cegval_labels_[elem.first] = lbl;
  }
}

template <class Prover_T>
void CegarValues<Prover_T>::cegar_abstract()
{
  prover_ts_ = conc_ts_;
  UnorderedTermMap prover_to_vals;
  ValueAbstractor va(prover_ts_, prover_to_vals);

  // now update with abstraction
  if (prover_ts_.is_functional()) {
    throw PonoException("Functional TS NYI in cegar_values");
  } else {
    Term init = prover_ts_.init();
    Term trans = prover_ts_.trans();
    static_cast<RelationalTransitionSystem &>(prover_ts_)
        .set_behavior(va.visit(init), va.visit(trans));
  }

  assert(super::bad_);
  super::bad_ = va.visit(super::bad_);
  cegval_bad_ = to_cegval_solver_.transfer_term(super::bad_, BOOL);

  // expecting to have had values to abstract
  assert(prover_to_vals.size());
  // copy over to the other solver
  for (const auto & elem : prover_to_vals) {
    // also make sure the abstract variables are frozen
    prover_ts_.assign_next(elem.first, elem.first);
    to_vals_[to_cegval_solver_.transfer_term(elem.first)] =
        to_cegval_solver_.transfer_term(elem.second);
  }
}

template <class Prover_T>
bool CegarValues<Prover_T>::cegar_refine()
{
  size_t cex_length = super::witness_length();

  // create bmc formula for abstract system
  Term bmcform = cegval_un_.at_time(cegval_ts_.init(), 0);
  for (size_t i = 0; i < cex_length; ++i) {
    bmcform = cegval_solver_->make_term(
        And, bmcform, cegval_un_.at_time(cegval_ts_.trans(), i));
  }
  bmcform = cegval_solver_->make_term(
      And, bmcform, cegval_un_.at_time(cegval_bad_, cex_length));

  cegval_solver_->push();
  cegval_solver_->assert_formula(bmcform);

  // TODO add lemmas to both cegval_ts_ and prover_ts_
  TermVec assumps;
  TermVec equalities;
  for (const auto & elem : to_vals_) {
    assert(cegval_ts_.is_curr_var(elem.first));
    assert(elem.second->is_value());
    Term lbl = cegval_labels_.at(elem.first);
    assumps.push_back(lbl);
    // since the variables are frozen, only need to add at time step 0
    Term eqval = cegval_solver_->make_term(Equal, elem.first, elem.second);
    equalities.push_back(eqval);
    Term eqval0 = cegval_un_.at_time(eqval, 0);
    Term imp = cegval_solver_->make_term(Implies, lbl, eqval0);
    cegval_solver_->assert_formula(imp);
  }
  assert(assumps.size() == equalities.size());

  Result r = cegval_solver_->check_sat_assuming(assumps);

  // do refinement if needed
  if (r.is_unsat()) {
    UnorderedTermSet core;
    cegval_solver_->get_unsat_core(core);
    for (size_t i = 0; i < assumps.size(); ++i) {
      if (core.find(assumps[i]) != core.end()) {
        logger.log(2, "CegarValues adding refinement axiom {}", equalities[i]);
        // need to refine both systems
        cegval_ts_.add_constraint(equalities[i]);
        Term transferred_eq =
            from_cegval_solver_.transfer_term(equalities[i], BOOL);
        // TODO this should be more modular
        //      can't assume super::abs_ts_ is the right one to constrain
        super::abs_ts_.add_constraint(transferred_eq);
      }
    }
  }
  cegval_solver_->pop();

  return r.is_unsat();
}

// TODO add other template classes
template class CegarValues<CegProphecyArrays<IC3IA>>;

}  // namespace pono
