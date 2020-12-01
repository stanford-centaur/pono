/*********************                                                  */
/*! \file ic3base.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan, Florian Lonsing
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Abstract base class implementation of IC3 parameterized by
**        the unit used in frames, pre-image computation, and inductive
**        and predecessor generalization techniques.
**
**/

#include "engines/ic3base.h"

#include <algorithm>

#include "assert.h"
#include "smt/available_solvers.h"
#include "utils/logger.h"

using namespace smt;
using namespace std;

namespace pono {

// helper functions

/** Less than comparison of the hash of two terms
 *  for use in sorting
 *  @param t0 the first term
 *  @param t1 the second term
 *  @return true iff t0's hash is less than t1's hash
 */
bool term_lt(const smt::Term & t0, const smt::Term & t1)
{
  return (t0->hash() < t1->hash());
}

/** IC3UnitHandler */
smt::Term IC3UnitHandler::smart_not(const Term & t) const
{
  Op op = t->get_op();
  if (op == Not) {
    TermVec children(t->begin(), t->end());
    assert(children.size() == 1);
    return children[0];
  } else {
    return solver_->make_term(Not, t);
  }
}

/** IC3Base */

IC3Base::IC3Base(Property & p, SolverEnum se, unique_ptr<IC3UnitHandler> && h)
    : super(p, se),
      handler_(std::move(h)),
      reducer_(create_solver(se)),
      solver_context_(0)
{
  initialize();
}

IC3Base::IC3Base(Property & p,
                 const SmtSolver & s,
                 unique_ptr<IC3UnitHandler> && h)
    : super(p, s),
      handler_(std::move(h)),
      reducer_(create_solver(s->get_solver_enum())),
      solver_context_(0)
{
  initialize();
}

IC3Base::IC3Base(const PonoOptions & opt,
                 Property & p,
                 SolverEnum se,
                 unique_ptr<IC3UnitHandler> && h)
    : super(opt, p, se),
      handler_(std::move(h)),
      reducer_(create_solver(se)),
      solver_context_(0)
{
  initialize();
}

IC3Base::IC3Base(const PonoOptions & opt,
                 Property & p,
                 const SmtSolver & s,
                 unique_ptr<IC3UnitHandler> && h)
    : super(opt, p, s),
      handler_(std::move(h)),
      reducer_(create_solver(s->get_solver_enum())),
      solver_context_(0)
{
  initialize();
}

void IC3Base::initialize()
{
  // TODO: fix initialize. Not sure it makes sense to call it here instead of in
  // prover constructor
  //       need to think about multiple levels of inheritance and where the
  //       responsibility for calling initialize belongs
  super::initialize();

  assert(solver_context_ == 0);  // expecting to be at base context level
  solver_true_ = solver_->make_term(true);

  frames_.clear();
  frame_labels_.clear();
  proof_goals_.clear();
  // first frame is always the initial states
  push_frame();
  // can't use constrain_frame for initial states because not guaranteed to be
  // an IC3Unit it's handled specially
  solver_->assert_formula(
      solver_->make_term(Implies, frame_labels_.at(0), ts_.init()));
  push_frame();

  // set semantics of TS labels
  Sort boolsort = solver_->make_sort(BOOL);
  assert(!init_label_);
  assert(!trans_label_);
  // frame 0 label is identical to init label
  init_label_ = frame_labels_[0];
  trans_label_ = solver_->make_symbol("__trans_label", boolsort);
  solver_->assert_formula(
      solver_->make_term(Implies, trans_label_, ts_.trans()));
}

ProverResult IC3Base::check_until(int k)
{
  check_ts();

  for (int i = 0; i <= k; ++i) {
    ProverResult r = step(i);
    if (r != ProverResult::UNKNOWN) {
      return r;
    }
  }

  return ProverResult::UNKNOWN;
}

bool IC3Base::witness(std::vector<smt::UnorderedTermMap> & out)
{
  throw PonoException("IC3 witness NYI");
}

// Protected Methods

bool IC3Base::intersects_bad()
{
  push_solver_context();
  // assert the last frame (conjunction over clauses)
  assert_frame_labels(reached_k_ + 1);
  // see if it intersects with bad
  solver_->assert_formula(bad_);
  Result r = solver_->check_sat();

  if (r.is_sat()) {
    // TODO: decide how important it is to push the whole bad for model based
    // IC3
    //       currently not doing that, would need to make intersects_bad virtual
    add_proof_goal(get_unit(), reached_k_ + 1, NULL);
  }

  pop_solver_context();

  assert(!r.is_unknown());
  return r.is_sat();
}

ProverResult IC3Base::step(int i)
{
  if (i <= reached_k_) {
    return ProverResult::UNKNOWN;
  }

  if (reached_k_ < 0) {
    return step_0();
  }

  // reached_k_ is the number of transitions that have been checked
  // at this point there are reached_k_ + 1 frames that don't
  // intersect bad, and reached_k_ + 2 frames overall
  assert(reached_k_ + 2 == frames_.size());
  logger.log(1, "Blocking phase at frame {}", i);
  // blocking phase
  while (intersects_bad()) {
    assert(has_proof_goals());
    if (!block_all()) {
      // counter-example
      return ProverResult::FALSE;
    }
  }

  logger.log(1, "Propagation phase at frame {}", i);
  // propagation phase
  push_frame();
  for (size_t j = 1; j < frames_.size() - 1; ++j) {
    if (propagate(j)) {
      assert(j + 1 < frames_.size());
      // save the invariant
      // which is the frame that just had all terms
      // from the previous frames propagated
      invar_ = get_frame(j + 1);
      return ProverResult::TRUE;
    }
  }

  ++reached_k_;

  return ProverResult::UNKNOWN;
}

ProverResult IC3Base::step_0()
{
  logger.log(1, "Checking if initial states satisfy property");
  assert(reached_k_ < 0);

  push_solver_context();
  solver_->assert_formula(init_label_);
  solver_->assert_formula(bad_);
  Result r = solver_->check_sat();
  if (r.is_sat()) {
    return ProverResult::FALSE;
  } else {
    assert(r.is_unsat());
    reached_k_ = 0;  // keep reached_k_ aligned with number of frames
  }
  pop_solver_context();
  return ProverResult::UNKNOWN;
}

bool IC3Base::get_predecessor(size_t i, const IC3Unit & c, IC3Unit & out_pred)
{
  assert(i > 0);
  assert(i < frames_.size());
  // expecting to be the polarity for proof goals, not frames
  // e.g. negated
  assert(c.negated);

  assert(solver_context_ == 0);
  push_solver_context();

  // F[i-1]
  assert_frame_labels(i - 1);
  // -c
  solver_->assert_formula(solver_->make_term(Not, c.term));
  // Trans
  assert_trans_label();
  // c'
  solver_->assert_formula(ts_.next(c.term));

  Result r = solver_->check_sat();
  if (r.is_sat()) {
    if (options_.ic3_pregen_) {
      out_pred = generalize_predecessor(i, c);
    } else {
      out_pred = get_unit();
    }
  } else {
    // TODO: consider automatically taking advantage
    //       of an unsat core. Took it out for now (was in MBIC3)
    //       because it needs to work for any IC3Unit
    //       Maybe IC3Unit needs to know how to generalize itself
    //         or at least how to make a conjunctive partition
    //         or it's possible they all can function approximately the same
    out_pred = c;
  }
  pop_solver_context();

  assert(!r.is_unknown());
  return r.is_sat();
}

// Helper methods

bool IC3Base::block_all()
{
  while (has_proof_goals()) {
    IC3Goal pg = get_next_proof_goal();
    // block can fail, which just means a
    // new proof goal will be added
    if (!block(pg) && !pg.idx) {
      // if a proof goal cannot be blocked at zero
      // then there's a counterexample
      return false;
    }
  }
  assert(!has_proof_goals());
  return true;
}

bool IC3Base::block(const IC3Goal & pg)
{
  const IC3Unit & c = pg.target;
  size_t i = pg.idx;

  logger.log(
      3, "Attempting to block proof goal <{}, {}>", c.term->to_string(), i);

  assert(i < frames_.size());
  assert(i >= 0);
  // TODO: assert c -> frames_[i]

  if (i == 0) {
    // can't block anymore -- this is a counterexample
    return false;
  }

  IC3Unit pred;  // populated by get_predecessor
  if (!get_predecessor(i, c, pred)) {
    assert(!pred.is_null());
    // can block this cube
    vector<IC3Unit> blocking_units;
    if (options_.ic3_indgen_) {
      blocking_units = inductive_generalization(i, pred);
    } else {
      IC3Unit bu = handler_->negate(pred);
      assert(bu.negated);
      blocking_units = { bu };
    }
    // pred is a subset of c
    logger.log(3, "Blocking term at frame {}: {}", i, c.term->to_string());
    if (options_.verbosity_ >= 3) {
      for (auto u : blocking_units) {
        logger.log(3, " with {}", u.term->to_string());
      }
    }

    // Most IC3 implementations will have only a single element in the vector
    // e.g. a single clause. But this is not guaranteed for all
    // for example, interpolant-based generalization for bit-vectors is not
    // always a single clause
    size_t min_idx = frames_.size();
    for (auto bu : blocking_units) {
      // TODO: fix name -- might not be a clause anymore
      // try to push
      size_t idx = find_highest_frame(i, bu);
      constrain_frame(idx, bu);
      if (idx < min_idx) {
        min_idx = idx;
      }
    }

    // we're limited by the minimum index that a conjunct could be pushed to
    if (min_idx + 1 < frames_.size()) {
      add_proof_goal(c, min_idx + 1, pg.next);
    }
    return true;
  } else {
    assert(!pred.is_null());
    add_proof_goal(pred, i - 1, make_shared<IC3Goal>(pg));
    return false;
  }
}

bool IC3Base::propagate(size_t i)
{
  assert(i + 1 < frames_.size());

  unordered_set<size_t> indices_to_remove;
  const vector<IC3Unit> & Fi = frames_.at(i);

  push_solver_context();
  assert_frame_labels(i);
  assert_trans_label();

  for (size_t j = 0; j < Fi.size(); ++j) {
    const Term & t = Fi.at(j).term;

    // Relative inductiveness check
    // Check F[i] /\ t /\ T /\ -t'
    // NOTE: asserting t is redundant because t \in F[i]
    push_solver_context();
    solver_->assert_formula(solver_->make_term(Not, ts_.next(t)));

    Result r = solver_->check_sat();
    assert(!r.is_unknown());
    if (r.is_unsat()) {
      // mark for removal
      indices_to_remove.insert(j);
      // add to next frame
      constrain_frame(i + 1, Fi.at(j));
    }

    pop_solver_context();
  }

  pop_solver_context();

  // keep invariant that terms are kept at highest frame
  // where they are known to hold
  // need to remove some terms from the current frame
  vector<IC3Unit> new_frame_i;
  new_frame_i.reserve(Fi.size() - indices_to_remove.size());
  for (size_t j = 0; j < Fi.size(); ++j) {
    if (indices_to_remove.find(j) == indices_to_remove.end()) {
      new_frame_i.push_back(Fi.at(j));
    }
  }

  frames_[i] = new_frame_i;

  return new_frame_i.empty();
}

void IC3Base::push_frame()
{
  assert(frame_labels_.size() == frames_.size());
  // pushes an empty frame
  frame_labels_.push_back(
      solver_->make_symbol("__frame_label_" + std::to_string(frames_.size()),
                           solver_->make_sort(BOOL)));
  frames_.push_back({});
}

void IC3Base::constrain_frame(size_t i, const IC3Unit & constraint)
{
  assert(i > 0);  // there's a special case for frame 0
  assert(i < frame_labels_.size());
  assert(frame_labels_.size() == frames_.size());
  solver_->assert_formula(
      solver_->make_term(Implies, frame_labels_.at(i), constraint.term));
  frames_.at(i).push_back(constraint);
}

void IC3Base::assert_frame_labels(size_t i) const
{
  // never expecting to assert a frame at base context
  assert(solver_context_ > 0);
  assert(frame_labels_.size() == frames_.size());
  Term assump;
  for (size_t j = 0; j < frame_labels_.size(); ++j) {
    assump = frame_labels_[j];
    if (j < i) {
      // optimization: disable the unused constraints
      // by asserting the negated label
      assump = solver_->make_term(Not, assump);
    }
    assert(assump);  // assert that it's non-null
    solver_->assert_formula(assump);
  }
}

Term IC3Base::get_frame(size_t i) const
{
  // TODO: decide if frames should hold IC3Units or terms
  //       need to special case initial state if using IC3Units
  if (i == 0) {
    // F[0] is always the initial states constraint
    return ts_.init();
  }

  Term res = solver_true_;
  for (size_t j = i; j < frames_.size(); ++j) {
    for (auto u : frames_[j]) {
      res = solver_->make_term(And, res, u.term);
    }
  }
  return res;
}

void IC3Base::assert_trans_label() const
{
  // shouldn't be a scenario where trans is asserted at base context
  // just because of how IC3 works
  assert(solver_context_ > 0);
  solver_->assert_formula(trans_label_);
}

bool IC3Base::has_proof_goals() const { return !proof_goals_.empty(); }

IC3Goal IC3Base::get_next_proof_goal()
{
  assert(has_proof_goals());
  IC3Goal pg = proof_goals_.back();
  proof_goals_.pop_back();
  return pg;
}

void IC3Base::add_proof_goal(const IC3Unit & c, size_t i, shared_ptr<IC3Goal> n)
{
  // IC3Unit aligned with frame so proof goal should be negated
  // e.g. for bit-level IC3, IC3Unit is a Clause and the proof
  // goal should be a Cube
  assert(c.negated);
  proof_goals_.push_back(IC3Goal(c, i, n));
}

bool IC3Base::intersects(const Term & A, const Term & B)
{
  push_solver_context();
  solver_->assert_formula(A);
  solver_->assert_formula(B);
  Result r = solver_->check_sat();
  pop_solver_context();
  return r.is_sat();
}

bool IC3Base::intersects_initial(const Term & t)
{
  return intersects(init_label_, t);
}

void IC3Base::fix_if_intersects_initial(TermVec & to_keep, const TermVec & rem)
{
  if (rem.size() != 0) {
    Term formula = solver_->make_term(And, init_label_, make_and(to_keep));
    reducer_.reduce_assump_unsatcore(formula,
                                     rem,
                                     to_keep,
                                     NULL,
                                     options_.ic3_gen_max_iter_,
                                     options_.random_seed_);
  }
}

size_t IC3Base::find_highest_frame(size_t i, const IC3Unit & u)
{
  assert(!u.negated);
  Term c = u.term;
  push_solver_context();
  solver_->assert_formula(c);
  solver_->assert_formula(solver_->make_term(Not, ts_.next(c)));
  assert_trans_label();

  Result r;
  size_t j;
  for (j = i; j + 1 < frames_.size(); ++j) {
    push_solver_context();
    assert_frame_labels(j);
    r = solver_->check_sat();
    pop_solver_context();
    if (r.is_sat()) {
      break;
    }
  }

  pop_solver_context();

  return j;
}

Term IC3Base::make_and(TermVec vec, SmtSolver slv) const
{
  if (!slv) {
    slv = solver_;
  }

  if (vec.size() == 0) {
    return slv->make_term(true);
  }

  // sort the conjuncts
  std::sort(vec.begin(), vec.end(), term_lt);
  Term res = vec[0];
  for (size_t i = 1; i < vec.size(); ++i) {
    res = slv->make_term(And, res, vec[i]);
  }
  return res;
}

void IC3Base::push_solver_context()
{
  solver_->push();
  solver_context_++;
}

void IC3Base::pop_solver_context()
{
  solver_->pop();
  solver_context_--;
}

}  // namespace pono
