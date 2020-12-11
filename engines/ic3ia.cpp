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

#include "smt-switch/cvc4_solver.h"
#include "smt-switch/cvc4_sort.h"
#include "smt-switch/cvc4_term.h"
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
      longest_cex_length_(0),
      cvc4_(create_solver(smt::SolverEnum::CVC4)),
      to_cvc4_(cvc4_),
      from_cvc4_(solver_)
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
      longest_cex_length_(0),
      cvc4_(create_solver(smt::SolverEnum::CVC4)),
      to_cvc4_(cvc4_),
      from_cvc4_(solver_)
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
      longest_cex_length_(0),
      cvc4_(create_solver(smt::SolverEnum::CVC4)),
      to_cvc4_(cvc4_),
      from_cvc4_(solver_)
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
      longest_cex_length_(0),
      cvc4_(create_solver(smt::SolverEnum::CVC4)),
      to_cvc4_(cvc4_),
      from_cvc4_(solver_)
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
  super::initialize();

  // add all the predicates from init and property to the abstraction
  // NOTE: abstract is called automatically in IC3Base initialize
  UnorderedTermSet preds;
  get_predicates(solver_, ts_->init(), preds, false);
  get_predicates(solver_, bad_, preds, false);
  for (auto p : preds) {
    add_predicate(p);
  }
  // more predicates will be added during refinement
  // these ones are just initial predicates

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
  // just need to set ts_ to the abstraction
  assert(abs_ts_.init());  // should be non-null
  assert(abs_ts_.trans());
  ts_ = &abs_ts_;
}

RefineResult IC3IA::refine()
{
  // recover the counterexample trace
  assert(intersects_initial(cex_pg_.target.term));
  TermVec cex({ cex_pg_.target.term });
  ProofGoal tmp = cex_pg_;
  while (tmp.next) {
    tmp = *(tmp.next);
    cex.push_back(tmp.target.term);
    assert(conc_ts_.only_curr(tmp.target.term));
  }
  assert(cex.size() > 1);

  UnorderedTermSet preds;

  // HACK added to experiment with CVC4 SyGuS for finding predicates
  if (options_.ic3ia_cvc4_pred_) {
    cvc4_find_preds(cex, preds);
  } else {
    // use interpolator to get predicates
    // remember -- need to transfer between solvers
    assert(interpolator_);
    Term t = make_and({ conc_ts_.init(), cex[0] });
    Term A = to_interpolator_.transfer_term(unroller_.at_time(t, 0), BOOL);

    TermVec B;
    B.reserve(cex.size() - 1);
    // add to B in reverse order so we can pop_back later
    for (int i = cex.size() - 1; i >= 0; --i) {
      register_symbol_mappings(
          i);  // make sure to_solver_ cache is populated with unrolled symbols
      t = cex[i];
      if (i + 1 < cex.size()) {
        t = solver_->make_term(And, t, conc_ts_.trans());
      }
      B.push_back(
          to_interpolator_.transfer_term(unroller_.at_time(t, i), BOOL));
    }

    // now get interpolants for each transition starting with the first
    bool all_sat = true;
    TermVec interpolants;
    while (B.size()) {
      // Note: have to pass the solver (defaults to solver_)
      Term fullB = make_and(B, interpolator_);
      Term I;
      Result r(smt::UNKNOWN, "unset result");
      try {
        r = interpolator_->get_interpolant(A, fullB, I);
      }
      catch (SmtException & e) {
        logger.log(3, e.what());
      }

      all_sat &= r.is_sat();
      if (r.is_unsat()) {
        Term untimedI = unroller_.untime(to_solver_.transfer_term(I, BOOL));
        logger.log(3, "got interpolant: {}", untimedI);
        interpolants.push_back(untimedI);
      }
      // move next cex time step to A
      // they were added to B in reverse order
      A = interpolator_->make_term(And, A, B.back());
      B.pop_back();
    }

    // record the length of this counterexample
    // important to set it here because it's used in register_symbol_mapping
    // to determine if state variables unrolled to a certain length
    // have already been cached in to_solver_
    longest_cex_length_ = cex.size();

    if (all_sat) {
      // this is a real counterexample, so the property is false
      return RefineResult::REFINE_NONE;
    } else if (!interpolants.size()) {
      logger.log(1,
                 "Interpolation failures...couldn't find any new predicates");
      return RefineResult::REFINE_FAIL;
    } else {
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

void IC3IA::register_symbol_mappings(size_t i)
{
  if (i < longest_cex_length_) {
    // these symbols should have already been handled
  }

  UnorderedTermMap & cache = to_solver_.get_cache();
  Term unrolled_sv;
  for (auto sv : ts_->statevars()) {
    unrolled_sv = unroller_.at_time(sv, i);
    cache[to_interpolator_.transfer_term(unrolled_sv)] = unrolled_sv;
  }
}

bool IC3IA::cvc4_find_preds(const TermVec & cex, UnorderedTermSet & out_preds)
{
  namespace cvc4a = ::CVC4::api;

  size_t cex_length = cex.size();
  assert(cex_length);

  // unroll the counterexample
  Term solver_formula = unroller_.at_time(abs_ts_.init(), 0);
  Term t;
  for (size_t i = 0; i < cex_length; ++i) {
    t = cex[i];
    if (i + 1 < cex_length) {
      t = solver_->make_term(And, t, abs_ts_.trans());
    }
    solver_formula =
        solver_->make_term(And, solver_formula, unroller_.at_time(t, i));
  }

  cvc4a::Term cvc4_formula = static_pointer_cast<CVC4Term>(
                                 to_cvc4_.transfer_term(solver_formula, BOOL))
                                 ->get_cvc4_term();
  vector<unordered_map<cvc4a::Term, cvc4a::Term, cvc4a::TermHashFunction>>
      cvc4_unroller;
  for (size_t i = 0; i < cex_length; ++i) {
    cvc4_unroller.push_back(
        unordered_map<cvc4a::Term, cvc4a::Term, cvc4a::TermHashFunction>());
    unordered_map<cvc4a::Term, cvc4a::Term, cvc4a::TermHashFunction> & m =
        cvc4_unroller.back();
    for (auto sv : conc_ts_.statevars()) {
      Term nv = conc_ts_.next(sv);
      Term abs_nv = ia_.abstract(nv);

      cvc4a::Term cvc4_sv = static_pointer_cast<CVC4Term>(sv)->get_cvc4_term();
      cvc4a::Term cvc4_nv = static_pointer_cast<CVC4Term>(nv)->get_cvc4_term();
      cvc4a::Term cvc4_abs_nv =
          static_pointer_cast<CVC4Term>(abs_nv)->get_cvc4_term();

      Term unrolled_sv = to_cvc4_.transfer_term(unroller_.at_time(sv, i));
      Term unrolled_nv = to_cvc4_.transfer_term(unroller_.at_time(nv, i));
      Term unrolled_abs_nv =
          to_cvc4_.transfer_term(unroller_.at_time(abs_nv, i));

      cvc4a::Term cvc4_unrolled_sv =
          static_pointer_cast<CVC4Term>(unrolled_sv)->get_cvc4_term();
      cvc4a::Term cvc4_unrolled_nv =
          static_pointer_cast<CVC4Term>(unrolled_nv)->get_cvc4_term();
      cvc4a::Term cvc4_unrolled_abs_nv =
          static_pointer_cast<CVC4Term>(unrolled_abs_nv)->get_cvc4_term();

      m[cvc4_sv] = cvc4_unrolled_sv;
      m[cvc4_nv] = cvc4_unrolled_nv;
      m[cvc4_abs_nv] = cvc4_unrolled_abs_nv;
    }
  }

  const cvc4a::Solver & cvc4_solver =
      static_pointer_cast<CVC4Solver>(cvc4_)->get_cvc4_solver();

  unordered_map<cvc4a::Term, cvc4a::Term, cvc4a::TermHashFunction> to_bound_var;
  unordered_map<cvc4a::Term, cvc4a::Term, cvc4a::TermHashFunction>
      from_bound_var;
  for (auto sv : conc_ts_.statevars()) {
    cvc4a::Term cvc4_sv = static_pointer_cast<CVC4Term>(sv)->get_cvc4_term();
    cvc4a::Term cvc4_bv =
        cvc4_solver.mkVar(cvc4_sv.getSort(), cvc4_sv.toString() + "_var");
    to_bound_var[cvc4_sv] = cvc4_bv;
    from_bound_var[cvc4_bv] = cvc4_sv;
  }

  // TODO finish this
  // create a synthFun, P, and add all the predicate unrollings to cvc4_formula
  // e.g. P(x^@0) <-> P(x@1) /\ P(x^@1) <-> P(x@2)
  throw PonoException("WIP");
}

}  // namespace pono
