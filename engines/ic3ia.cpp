/*********************                                                  */
/*! \file ic3ia.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
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
**/

#include "engines/ic3ia.h"

#include "utils/logger.h"
#include "utils/term_analysis.h"

using namespace smt;
using namespace std;

namespace pono {

IC3IA::IC3IA(Property & p, SolverEnum se)
    : super(p, se), abs_ts_(ts_.solver()), ia_(ts_, abs_ts_)
{
  initialize();
}

IC3IA::IC3IA(Property & p, const SmtSolver & slv)
    : super(p, slv), abs_ts_(ts_.solver()), ia_(ts_, abs_ts_)
{
  initialize();
}

IC3IA::IC3IA(const PonoOptions & opt, Property & p, const SolverEnum se)
    : super(opt, p, se), abs_ts_(ts_.solver()), ia_(ts_, abs_ts_)
{
  initialize();
}

IC3IA::IC3IA(const PonoOptions & opt, Property & p, const SmtSolver & slv)
    : super(opt, p, slv), abs_ts_(ts_.solver()), ia_(ts_, abs_ts_)
{
  initialize();
}

IC3IA::~IC3IA() {}

// protected methods
void IC3IA::initialize()
{
  boolsort_ = solver_->make_sort(BOOL);

  if (options_.ic3_indgen_mode_ == 2) {
    // TODO: clean this up
    throw PonoException(
        "Interpolant-based generalization not supported with IC3IA -- doesn't "
        "really make sense");
  }

  // set up the interpolator for refinement
  assert(options_.ic3_indgen_mode_ != 2);
  assert(!interpolator_);
  assert(!to_interpolator_);
  assert(!to_solver_);
  initialize_interpolator();

  interp_ts_ = RelationalTransitionSystem(ts_, *to_interpolator_);
  interp_unroller_ = make_unique<Unroller>(interp_ts_, interpolator_);
}

bool IC3IA::block_all()
{
  while (has_proof_goals()) {
    ProofGoal pg = get_next_proof_goal();
    // block can fail, which just means a
    // new proof goal will be added
    if (!block(pg) && !pg.idx) {
      // if a proof goal cannot be blocked at zero
      // then there's a (possibly abstract) counterexample
      ProverResult refined = refine(pg);
      if (refined == ProverResult::FALSE) {
        // got a real counterexample trace
        return false;
      } else if (refined == ProverResult::UNKNOWN) {
        // TODO: feed this through so it returns unknown
        throw PonoException("Could not refine");
      }
    }
  }
  assert(!has_proof_goals());
  return true;
}

bool IC3IA::intersects_bad()
{
  // TODO: if we refactor ModelBasedIC3 to use get_conjunction_from_model
  //       then we shouldn't have to override this at all

  assert(reached_k_ + 2 == frames_.size());
  assert(frame_labels_.size() == frames_.size());

  push_solver_context();

  // check if last frame intersects with bad
  assert_frame(reached_k_ + 1);
  solver_->assert_formula(bad_);
  Result r = solver_->check_sat();

  if (r.is_sat()) {
    add_proof_goal(get_conjunction_from_model(), reached_k_ + 1, NULL);
  }

  pop_solver_context();

  return r.is_sat();
}

bool IC3IA::get_predecessor(size_t i,
                            const Conjunction & c,
                            Conjunction & out_pred)
{
  assert(i > 0);
  assert(i < frames_.size());

  assert(solver_context_ == 0);
  push_solver_context();

  // F[i-1]
  assert_frame(i - 1);
  // -c
  solver_->assert_formula(solver_->make_term(Not, c.term_));
  // Trans
  solver_->assert_formula(trans_label_);
  // c'
  solver_->assert_formula(abs_ts_.next(c.term_));

  Result r = solver_->check_sat();
  if (r.is_sat()) {
    // TODO: add generalize option
    // TODO: refactor mbic3 to use get_conjunction_from_model and we won't need
    // to override it
    out_pred = get_conjunction_from_model();
    pop_solver_context();
  } else {
    pop_solver_context();
    assert(solver_context_ == 0);
    TermVec assump, red_assump, rem_assump;
    for (auto a : c.conjuncts_) {
      assump.push_back(abs_ts_.next(a));
    }

    Term formula = make_and(TermVec{
        get_frame(i - 1), solver_->make_term(Not, c.term_), trans_label_ });

    // filter using unsatcore
    reduce_assump_unsatcore(formula, assump, red_assump, &rem_assump);
    // get current version of red_assump
    TermVec cur_red_assump, cur_rem_assump;
    for (auto a : red_assump) {
      cur_red_assump.push_back(abs_ts_.curr(a));
    }
    for (auto a : rem_assump) {
      cur_rem_assump.push_back(abs_ts_.curr(a));
    }

    fix_if_intersects_initial(cur_red_assump, cur_rem_assump);
    assert(cur_red_assump.size() > 0);
    out_pred = Conjunction(solver_, cur_red_assump);
  }

  assert(!r.is_unknown());
  return r.is_sat();
}

Conjunction IC3IA::get_conjunction_from_model()
{
  TermVec conjuncts;
  conjuncts.reserve(pred_statevars_.size());
  Term val;
  for (auto p : pred_statevars_) {
    if ((val = solver_->get_value(p)) == true_) {
      conjuncts.push_back(p);
    } else {
      assert(val == false_);
      conjuncts.push_back(solver_->make_term(Not, p));
    }
  }
  return Conjunction(solver_, conjuncts);
}

void IC3IA::set_labels()
{
  assert(solver_context_ == 0);  // expecting to be at base context level
  // set semantics of labels
  Sort boolsort = solver_->make_sort(BOOL);
  if (!init_label_) {
    // frame 0 label is identical to init label
    init_label_ = frame_labels_[0];
  }
  if (!trans_label_) {
    trans_label_ = solver_->make_symbol("__trans_label", boolsort);
    // Note: Using iff instead of implies so that we can negate it
    // for example in relational preimage generalization
    solver_->assert_formula(
        solver_->make_term(Iff, trans_label_, abs_ts_.trans()));

    // add boolean state variables as predicates
    // NOTE: this must be done before adding predicates
    // because adding predicates also creates a new state var
    // for each predicate and we don't want to dobule count those
    for (auto sv : abs_ts_.statevars()) {
      if (boolsort_ == sv->get_sort()) {
        pred_statevars_.push_back(sv);
      }
    }

    // add all the predicates from init and property
    UnorderedTermSet preds;
    get_predicates(ts_.init(), boolsort_, preds, false);
    get_predicates(bad_, boolsort_, preds, false);
    for (auto p : preds) {
      add_predicate(p);
    }
  }
}

bool IC3IA::only_curr(Term & t) { return abs_ts_.only_curr(t); }

Term IC3IA::next(const Term & t) const { return abs_ts_.next(t); }

void IC3IA::set_invar(size_t i)
{
  Term Fi = get_frame(i);
  // need to translate to actual predicates
  invar_ = solver_->substitute(Fi, predvar2pred_);
}

void IC3IA::add_predicate(const Term & pred)
{
  assert(abs_ts_.only_curr(pred));
  // add predicate to abstraction and get the new constraint
  Term predabs_rel = ia_.add_predicate(pred);
  // refine the transition relation incrementally
  // by adding a new constraint
  assert(!solver_context_);  // should be at context 0
  solver_->assert_formula(
      solver_->make_term(Implies, trans_label_, predabs_rel));

  // create a new state variable for the IC3 procedure
  // used to represent the predicate
  // NOTE: don't need to actually add to the TS
  //       kept separate for predicate abstraction
  // TODO: add infrastructure for getting a fresh symbol -- might fail if name
  // already exists
  string name = "__pred" + std::to_string(pred_statevars_.size());
  Term pred_state = abs_ts_.make_statevar(name, boolsort_);
  pred_statevars_.push_back(pred_state);
  predvar2pred_[pred_state] = pred;

  // associate the fresh state vars with predicate applied to current / next
  // states
  // NOTE: needs to be Iff, not just Implies, because needs to enforce false
  // case also
  solver_->assert_formula(solver_->make_term(Iff, pred_state, pred));
  // NOTE: this predicate is over the original next state vars (not abstracted)
  solver_->assert_formula(
      solver_->make_term(Iff, abs_ts_.next(pred_state), ts_.next(pred)));
}

ProverResult IC3IA::refine(ProofGoal pg)
{
  // recover the counterexample trace
  assert(intersects_initial(pg.conj.term_));
  TermVec cex({ pg.conj.term_ });
  ProofGoal tmp = pg;
  while (tmp.next) {
    tmp = *(tmp.next);
    cex.push_back(tmp.conj.term_);
  }
  assert(cex.size() > 1);

  // use interpolator to get predicates
  // remember -- need to transfer between solvers
  assert(interpolator_);
  Term t = make_and({ ts_.init(), ts_.trans(), cex[0] });
  Term A =
      interp_unroller_->at_time(to_interpolator_->transfer_term(t, BOOL), 0);

  TermVec B;
  Term interp_trans = to_interpolator_->transfer_term(ts_.trans(), BOOL);
  B.reserve(cex.size() - 1);
  // add to B in reverse order so we can pop_back later
  for (int i = cex.size() - 1; i >= 0; --i) {
    // replace indicator variables with actual predicates
    Term cex_preds = solver_->substitute(cex[i], predvar2pred_);
    t = to_interpolator_->transfer_term(cex_preds, BOOL);
    if (i + 1 < cex.size()) {
      t = interpolator_->make_term(And, t, interp_trans);
    }
    B.push_back(interp_unroller_->at_time(t, i));
  }

  // now get interpolants for each transition starting with the first
  bool all_sat = true;
  TermVec interpolants;
  while (B.size()) {
    // Note: have to pass the solver (defaults to solver_)
    Term fullB = make_and(B, interpolator_);
    Term I;
    try {
      Result r = interpolator_->get_interpolant(A, fullB, I);
      all_sat &= r.is_sat();
      if (r.is_unsat()) {
        Term untimedI = interp_unroller_->untime(I);
        logger.log(3, "got interpolant: {}", untimedI);
        interpolants.push_back(to_solver_->transfer_term(untimedI, BOOL));
      }
    }
    catch (SmtException & e) {
      logger.log(3, e.what());
    }
    // move next cex time step to A
    // they were added to B in reverse order
    A = interpolator_->make_term(And, A, B.back());
    B.pop_back();
  }

  if (all_sat) {
    // this is a real counterexample, so the property is false
    return ProverResult::FALSE;
  } else if (!interpolants.size()) {
    logger.log(1, "Interpolation failures...couldn't find any new predicates");
    return ProverResult::UNKNOWN;
  } else {
    UnorderedTermSet preds;
    for (auto I : interpolants) {
      assert(ts_.only_curr(I));
      get_predicates(I, boolsort_, preds);
    }

    if (!preds.size()) {
      logger.log(1, "No new predicates found...");
      return ProverResult::UNKNOWN;
    }

    // add all the new predicates
    for (auto p : preds) {
      add_predicate(p);
    }

    // able to refine the system to rule out this abstract counterexample
    return ProverResult::TRUE;
  }
}

}  // namespace pono
