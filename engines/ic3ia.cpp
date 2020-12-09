/*********************                                                  */
/*! \file ic3ia.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief IC3 via Implicit Predicate Abstraction (IC3IA) implementation
**        based on
**
**        IC3 Modulo Theories via Implicit Predicate Abstraction
**            -- Alessandro Cimatti, Alberto Griggio,
**               Sergio Mover, Stefano Tonetta
**
**        and the open source implementation:
**
**        https://es-static.fbk.eu/people/griggio/ic3ia/index.html
**
**  within Pono, we are building on the bit-level IC3 instead of directly
**  on IC3Base, because a lot of the functionality is the same
**  In particular, we don't need to override either of the generalization
**  functions. Instead focusing on abstract/refine.
**
**/

#include "engines/ic3ia.h"

#include <random>

#include "smt/available_solvers.h"
#include "utils/logger.h"
#include "utils/term_analysis.h"

using namespace smt;
using namespace std;

namespace pono {

IC3IA::IC3IA(Property & p, SolverEnum se, SolverEnum itp_se)
    : super(p, se),
      conc_ts_(property_.transition_system()),
      abs_ts_(solver_),
      ia_(conc_ts_, abs_ts_, unroller_),
      interpolator_(create_interpolating_solver(itp_se)),
      to_interpolator_(interpolator_),
      to_solver_(solver_),
      interp_ts_(conc_ts_, to_interpolator_),
      interp_unroller_(interp_ts_, interpolator_)
{
}

IC3IA::IC3IA(Property & p, const SmtSolver & s, SolverEnum itp_se)
    : super(p, s),
      conc_ts_(property_.transition_system()),
      abs_ts_(solver_),
      ia_(conc_ts_, abs_ts_, unroller_),
      interpolator_(create_interpolating_solver(itp_se)),
      to_interpolator_(interpolator_),
      to_solver_(solver_),
      interp_ts_(conc_ts_, to_interpolator_),
      interp_unroller_(interp_ts_, interpolator_)
{
}

IC3IA::IC3IA(const PonoOptions & opt,
             Property & p,
             SolverEnum se,
             SolverEnum itp_se)
    : super(opt, p, se),
      conc_ts_(property_.transition_system()),
      abs_ts_(solver_),
      ia_(conc_ts_, abs_ts_, unroller_),
      interpolator_(create_interpolating_solver(itp_se)),
      to_interpolator_(interpolator_),
      to_solver_(solver_),
      interp_ts_(conc_ts_, to_interpolator_),
      interp_unroller_(interp_ts_, interpolator_)
{
}

IC3IA::IC3IA(const PonoOptions & opt,
             Property & p,
             const SmtSolver & s,
             SolverEnum itp_se)
    : super(opt, p, s),
      conc_ts_(property_.transition_system()),
      abs_ts_(solver_),
      ia_(conc_ts_, abs_ts_, unroller_),
      interpolator_(create_interpolating_solver(itp_se)),
      to_interpolator_(interpolator_),
      to_solver_(solver_),
      interp_ts_(conc_ts_, to_interpolator_),
      interp_unroller_(interp_ts_, interpolator_)
{
}

// pure virtual method implementations

IC3Formula IC3IA::get_ic3_formula(TermVec * inputs, TermVec * nexts) const
{
  const TermVec & preds = ia_.predicates();
  TermVec conjuncts;
  conjuncts.reserve(preds.size());
  for (auto p : preds) {
    if (solver_->get_value(p) == solver_true_) {
      conjuncts.push_back(p);
    } else {
      conjuncts.push_back(solver_->make_term(Not, p));
    }

    if (nexts) {
      Term next_p = ts_->next(p);
      if (solver_->get_value(next_p) == solver_true_) {
        nexts->push_back(next_p);
      } else {
        nexts->push_back(solver_->make_term(Not, next_p));
      }
    }
  }

  if (inputs) {
    for (auto iv : ts_->inputvars()) {
      inputs->push_back(solver_->make_term(Equal, iv, solver_->get_value(iv)));
    }
  }

  return ic3_formula_conjunction(conjuncts);
}

bool IC3IA::ic3_formula_check_valid(const IC3Formula & u) const
{
  Sort boolsort = solver_->make_sort(BOOL);
  // check that children are literals
  Term pred;
  Op op;
  for (auto c : u.children) {
    if (c->get_sort() != boolsort) {
      logger.log(3, "ERROR IC3IA IC3Formula contains non-boolean atom: {}", c);
      return false;
    }

    pred = c;
    op = pred->get_op();
    if (op == Not || op == BVNot) {
      pred = *(c->begin());
      assert(pred);
    }

    // expecting either a boolean variable or a predicate
    if (!pred->is_symbolic_const() && predset_.find(pred) == predset_.end()) {
      logger.log(3, "ERROR IC3IA IC3Formula contains unknown atom: {}", pred);
      return false;
    }
  }

  // got through all checks without failing
  return true;
}

void IC3IA::check_ts() const
{
  // basically a No-Op
  // no restrictions except that interpolants must be supported
  // instead of checking explicitly, just let the interpolator throw an
  // exception better than maintaining in two places
}

void IC3IA::initialize()
{
  // set ts_ to the abstraction BEFORE initializing base classes
  ts_ = &abs_ts_;

  super::initialize();

  // populate cache for existing terms in solver_
  UnorderedTermMap & cache = to_solver_.get_cache();
  Term ns;
  for (auto s : ts_->statevars()) {
    // common variables are next states, unless used for refinement in IC3IA
    // then will refer to current state variables after untiming
    // need to cache both
    cache[to_interpolator_.transfer_term(s)] = s;
    ns = ts_->next(s);
    cache[to_interpolator_.transfer_term(ns)] = ns;
  }

  // need to add uninterpreted functions as well
  // first need to find them all
  // NOTE need to use get_free_symbols NOT get_free_symbolic_consts
  // because the latter ignores uninterpreted functions
  UnorderedTermSet free_symbols;
  get_free_symbols(ts_->init(), free_symbols);
  get_free_symbols(ts_->trans(), free_symbols);
  get_free_symbols(bad_, free_symbols);

  for (auto s : free_symbols) {
    assert(s->is_symbol());
    if (s->is_symbolic_const()) {
      // ignore constants
      continue;
    }
    cache[to_interpolator_.transfer_term(s)] = s;
  }

  // TODO fix generalize_predecessor for ic3ia
  //      might need to override it
  //      behaves a bit differently with both concrete and abstract next state
  //      vars
  if (options_.ic3_pregen_) {
    logger.log(1,
               "WARNING automatically disabling predecessor generalization -- "
               "not supported in IC3IA yet.");
    options_.ic3_pregen_ = false;
  }
}

void IC3IA::abstract()
{
  // main abstraction already done in constructor of ia_

  // add all the predicates from init and property
  UnorderedTermSet preds;
  get_predicates(solver_, ts_->init(), preds, false);
  get_predicates(solver_, bad_, preds, false);
  for (auto p : preds) {
    add_predicate(p);
  }
  // more predicates will be added during refinement
  // these ones are just initial predicates
}

RefineResult IC3IA::refine()
{
  // recover the counterexample trace
  assert(intersects_initial(cex_pg_.target.term));
  TermVec cex({ cex_pg_.target.term });
  IC3Goal tmp = cex_pg_;
  while (tmp.next) {
    tmp = *(tmp.next);
    cex.push_back(tmp.target.term);
    assert(conc_ts_.only_curr(tmp.target.term));
  }
  assert(cex.size() > 1);

  // use interpolator to get predicates
  // remember -- need to transfer between solvers
  assert(interpolator_);
  Term t = make_and({ conc_ts_.init(), cex[0] });
  Term A = interp_unroller_.at_time(to_interpolator_.transfer_term(t, BOOL), 0);

  TermVec B;
  Term interp_trans = to_interpolator_.transfer_term(conc_ts_.trans(), BOOL);
  B.reserve(cex.size() - 1);
  // add to B in reverse order so we can pop_back later
  for (int i = cex.size() - 1; i >= 0; --i) {
    // replace indicator variables with actual predicates
    t = to_interpolator_.transfer_term(cex[i], BOOL);
    if (i + 1 < cex.size()) {
      t = interpolator_->make_term(And, t, interp_trans);
    }
    B.push_back(interp_unroller_.at_time(t, i));
  }

  // now get interpolants for each transition starting with the first
  bool all_sat = true;
  TermVec interpolants;
  while (B.size()) {
    // Note: have to pass the solver (defaults to solver_)
    Term fullB = make_and(B, interpolator_);
    Term I;
    Result r(smt::UNKNOWN, "not set");
    try {
      r = interpolator_->get_interpolant(A, fullB, I);
    }
    catch (SmtException & e) {
      logger.log(3, e.what());
    }

    all_sat &= r.is_sat();
    if (r.is_unsat()) {
      Term untimedI = interp_unroller_.untime(I);
      logger.log(3, "got interpolant: {}", untimedI);
      interpolants.push_back(to_solver_.transfer_term(untimedI, BOOL));
    }
    // move next cex time step to A
    // they were added to B in reverse order
    A = interpolator_->make_term(And, A, B.back());
    B.pop_back();
  }

  if (all_sat) {
    // this is a real counterexample, so the property is false
    return RefineResult::REFINE_NONE;
  } else if (!interpolants.size()) {
    logger.log(1, "Interpolation failures...couldn't find any new predicates");
    return RefineResult::REFINE_FAIL;
  } else {
    UnorderedTermSet preds;
    for (auto I : interpolants) {
      assert(conc_ts_.only_curr(I));
      get_predicates(solver_, I, preds);
    }

    // reduce new predicates
    TermVec preds_vec(preds.begin(), preds.end());
    if (options_.random_seed_ > 0) {
      shuffle(preds_vec.begin(),
              preds_vec.end(),
              default_random_engine(options_.random_seed_));
    }

    TermVec red_preds;
    if (ia_.reduce_predicates(cex, preds_vec, red_preds)) {
      // reduction successful
      preds.clear();
      preds.insert(red_preds.begin(), red_preds.end());
    }

    // add all the new predicates
    bool found_new_preds = false;
    for (auto p : preds) {
      found_new_preds |= add_predicate(p);
    }

    if (!found_new_preds) {
      logger.log(1, "No new predicates found...");
      return RefineResult::REFINE_FAIL;
    }

    // able to refine the system to rule out this abstract counterexample
    return RefineResult::REFINE_SUCCESS;
  }
}

bool IC3IA::add_predicate(const Term & pred)
{
  if (predset_.find(pred) != predset_.end()) {
    // don't allow re-adding the same predicate
    return false;
  }

  assert(ts_->only_curr(pred));
  logger.log(2, "adding predicate {}", pred);
  predset_.insert(pred);
  // add predicate to abstraction and get the new constraint
  Term predabs_rel = ia_.add_predicate(pred);
  // refine the transition relation incrementally
  // by adding a new constraint
  assert(!solver_context_);  // should be at context 0
  solver_->assert_formula(
      solver_->make_term(Implies, trans_label_, predabs_rel));
  return true;
}

}  // namespace pono
