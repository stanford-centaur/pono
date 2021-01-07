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

IC3IA::IC3IA(const Property & p,
             const TransitionSystem & ts,
             const SmtSolver & s,
             PonoOptions opt)
    : super(p, RelationalTransitionSystem(s), s, opt),
      conc_ts_(ts, to_prover_solver_),
      ia_(conc_ts_, ts_, unroller_),
      // only mathsat interpolator supported
      interpolator_(create_interpolating_solver_for(
          SolverEnum::MSAT_INTERPOLATOR, Engine::IC3IA_ENGINE)),
      to_interpolator_(interpolator_),
      to_solver_(solver_),
      longest_cex_length_(0)
{
  engine_ = Engine::IC3IA_ENGINE;
}

// pure virtual method implementations

IC3Formula IC3IA::get_model_ic3formula(TermVec * out_inputs,
                                       TermVec * out_nexts) const
{
  const TermVec & preds = ia_.predicates();
  TermVec conjuncts;
  conjuncts.reserve(preds.size());
  for (const auto &p : preds) {
    if (solver_->get_value(p) == solver_true_) {
      conjuncts.push_back(p);
    } else {
      conjuncts.push_back(solver_->make_term(Not, p));
    }

    if (out_nexts) {
      Term next_p = ts_.next(p);
      if (solver_->get_value(next_p) == solver_true_) {
        out_nexts->push_back(next_p);
      } else {
        out_nexts->push_back(solver_->make_term(Not, next_p));
      }
    }
  }

  if (out_inputs) {
    for (const auto &iv : ts_.inputvars()) {
      out_inputs->push_back(
          solver_->make_term(Equal, iv, solver_->get_value(iv)));
    }
  }

  return ic3formula_conjunction(conjuncts);
}

bool IC3IA::ic3formula_check_valid(const IC3Formula & u) const
{
  const Sort &boolsort = solver_->make_sort(BOOL);
  // check that children are literals
  Term pred;
  Op op;
  for (const auto &c : u.children) {
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
    if (predset_.find(pred) == predset_.end()) {
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
  super::initialize();

  Sort boolsort = solver_->make_sort(BOOL);
  // add predicates automatically added by ia_
  // to our predset_
  // needed to prevent adding duplicate predicates later
  for (const auto & p : ia_.predicates()) {
    // ia_ automatically includes all boolean variables
    assert(p->is_symbolic_const() && p->get_sort() == boolsort);
    predset_.insert(p);
  }

  // add all the predicates from init and property to the abstraction
  // NOTE: abstract is called automatically in IC3Base initialize
  UnorderedTermSet preds;
  get_predicates(solver_, conc_ts_.init(), preds, false, false, true);
  size_t num_init_preds = preds.size();
  get_predicates(solver_, bad_, preds, false, false, true);
  size_t num_prop_preds = preds.size() - num_init_preds;
  for (const auto &p : preds) {
    add_predicate(p);
  }
  logger.log(1, "Number predicates found in init: {}", num_init_preds);
  logger.log(1, "Number predicates found in prop: {}", num_prop_preds);
  logger.log(1, "Total number of initial predicates: {}", preds.size());
  // more predicates will be added during refinement
  // these ones are just initial predicates

  // populate cache for existing terms in solver_
  UnorderedTermMap & cache = to_solver_.get_cache();
  Term ns;
  for (auto const&s : ts_.statevars()) {
    // common variables are next states, unless used for refinement in IC3IA
    // then will refer to current state variables after untiming
    // need to cache both
    cache[to_interpolator_.transfer_term(s)] = s;
    ns = ts_.next(s);
    cache[to_interpolator_.transfer_term(ns)] = ns;
  }

  // need to add uninterpreted functions as well
  // first need to find them all
  // NOTE need to use get_free_symbols NOT get_free_symbolic_consts
  // because the latter ignores uninterpreted functions
  UnorderedTermSet free_symbols;
  get_free_symbols(ts_.init(), free_symbols);
  get_free_symbols(ts_.trans(), free_symbols);
  get_free_symbols(bad_, free_symbols);

  for (auto const&s : free_symbols) {
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
  ia_.do_abstraction();

  assert(ts_.init());  // should be non-null
  assert(ts_.trans());
}

RefineResult IC3IA::refine()
{
  // recover the counterexample trace
  assert(check_intersects_initial(cex_pg_->target.term));
  TermVec cex;
  const ProofGoal * tmp = cex_pg_;
  while (tmp) {
    cex.push_back(tmp->target.term);
    assert(ts_.only_curr(tmp->target.term));
    tmp = tmp->next;
  }

  if (cex.size() == 1) {
    // if there are no transitions, then this is a concrete CEX
    return REFINE_NONE;
  }

  size_t cex_length = cex.size();

  // use interpolator to get predicates
  // remember -- need to transfer between solvers
  assert(interpolator_);

  TermVec formulae;
  for (size_t i = 0; i < cex_length; ++i) {
    // make sure to_solver_ cache is populated with unrolled symbols
    register_symbol_mappings(i);

    Term t = unroller_.at_time(cex[i], i);
    if (i + 1 < cex_length) {
      t = solver_->make_term(And, t, unroller_.at_time(conc_ts_.trans(), i));
    }
    formulae.push_back(to_interpolator_.transfer_term(t, BOOL));
  }

  TermVec out_interpolants;
  Result r =
      interpolator_->get_sequence_interpolants(formulae, out_interpolants);

  if (r.is_sat()) {
    // this is a real counterexample, so the property is false
    return RefineResult::REFINE_NONE;
  }

  // record the length of this counterexample
  // important to set it here because it's used in register_symbol_mapping
  // to determine if state variables unrolled to a certain length
  // have already been cached in to_solver_
  longest_cex_length_ = cex.size();

  UnorderedTermSet preds;
  for (auto const&I : out_interpolants) {
    if (!I) {
      assert(
          r.is_unknown());  // should only have null terms if got unknown result
      continue;
    }

    Term solver_I = unroller_.untime(to_solver_.transfer_term(I, BOOL));
    assert(conc_ts_.only_curr(solver_I));
    logger.log(3, "got interpolant: {}", solver_I);
    get_predicates(solver_, solver_I, preds, false, false, true);
  }

  // new predicates
  TermVec fresh_preds;
  for (auto const&p : preds) {
    if (predset_.find(p) == predset_.end()) {
      // unseen predicate
      fresh_preds.push_back(p);
    }
  }

  if (!fresh_preds.size()) {
    logger.log(1, "IC3IA: refinement failed couldn't find any new predicates");
    return RefineResult::REFINE_FAIL;
  }

  if (options_.random_seed_ > 0) {
    shuffle(fresh_preds.begin(),
            fresh_preds.end(),
            default_random_engine(options_.random_seed_));
  }

  // reduce new predicates
  TermVec red_preds;
  if (ia_.reduce_predicates(cex, fresh_preds, red_preds)) {
    // reduction successful
    logger.log(2,
               "reduce predicates successful {}/{}",
               red_preds.size(),
               fresh_preds.size());
    if (red_preds.size() < fresh_preds.size()) {
      fresh_preds.clear();
      fresh_preds.insert(fresh_preds.end(), red_preds.begin(), red_preds.end());
    }
  } else {
    // should only fail if removed all predicates
    // this can happen when there are uninterpreted functions
    // the unrolling can force incompatible UF interpretations
    // but IC3 (which doesn't unroll) still needs the predicates
    // in this case, just use all the fresh predicates
    assert(red_preds.size() == 0);
    logger.log(2, "reduce predicates FAILED");
  }

  // add all the new predicates
  for (auto const & p : fresh_preds) {
    bool new_pred = add_predicate(p);
    // expect all predicates to be new (e.g. unseen)
    // they were already filtered above
    assert(new_pred);
  }

  logger.log(1, "{} new predicates added by refinement", fresh_preds.size());

  // able to refine the system to rule out this abstract counterexample
  return RefineResult::REFINE_SUCCESS;
}

bool IC3IA::add_predicate(const Term & pred)
{
  if (predset_.find(pred) != predset_.end()) {
    // don't allow re-adding the same predicate
    return false;
  }

  assert(ts_.only_curr(pred));
  logger.log(2, "adding predicate {}", pred);
  predset_.insert(pred);
  // add predicate to abstraction and get the new constraint
  Term predabs_rel = ia_.add_predicate(pred);
  // refine the transition relation incrementally
  // by adding a new constraint
  assert(!solver_context_);  // should be at context 0
  solver_->assert_formula(
      solver_->make_term(Implies, trans_label_, predabs_rel));

  assert(predset_.size() == ia_.predicates().size());

  return true;
}

void IC3IA::register_symbol_mappings(size_t i)
{
  if (i < longest_cex_length_) {
    // these symbols should have already been handled
  }

  UnorderedTermMap & cache = to_solver_.get_cache();
  Term unrolled_sv;
  for (const auto &sv : ts_.statevars()) {
    unrolled_sv = unroller_.at_time(sv, i);
    cache[to_interpolator_.transfer_term(unrolled_sv)] = unrolled_sv;
  }
}

}  // namespace pono
