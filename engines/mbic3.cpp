/*********************                                         */
/*! \file mbic3.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan, Florian Lonsing
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Simple implementation of IC3 operating on a functional
**        transition system (exploiting this structure for
**        predecessor computation) and uses models.
**/
#include "engines/mbic3.h"
#include "utils/logger.h"

using namespace smt;
using namespace std;

namespace pono {

ModelBasedIC3::ModelBasedIC3(const Property & p, smt::SmtSolver & slv)
    : super(p, slv),
      true_(solver_->make_term(true)),
      false_(solver_->make_term(false))
{
  initialize();
}

ModelBasedIC3::ModelBasedIC3(const PonoOptions & opt,
                             const Property p,
                             smt::SmtSolver & slv)
    : super(opt, p, slv),
      true_(solver_->make_term(true)),
      false_(solver_->make_term(false))
{
  initialize();
}

ModelBasedIC3::~ModelBasedIC3() {}

void ModelBasedIC3::initialize()
{
  super::initialize();
  frames_.clear();
  proof_goals_.clear();
  // first frame is always the initial states
  frames_.push_back({ ts_.init() });
  push_frame();
}

ProverResult ModelBasedIC3::check_until(int k)
{
  for (int i = 0; i <= k; ++i) {
    ProverResult r = step(i);
    if (r != ProverResult::UNKNOWN) {
      return r;
    }
  }

  return ProverResult::UNKNOWN;
}

ProverResult ModelBasedIC3::step(int i)
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
    assert(!proof_goals_.empty());
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
      return ProverResult::TRUE;
    }
  }

  ++reached_k_;

  return ProverResult::UNKNOWN;
}

ProverResult ModelBasedIC3::step_0()
{
  logger.log(1, "Checking if initial states satisfy property");
  assert(reached_k_ < 0);

  solver_->push();
  solver_->assert_formula(ts_.init());
  solver_->assert_formula(bad_);
  Result r = solver_->check_sat();
  if (r.is_sat()) {
    return ProverResult::FALSE;
  } else {
    assert(r.is_unsat());
    reached_k_ = 0;  // keep reached_k_ aligned with number of frames
  }
  solver_->pop();
  return ProverResult::UNKNOWN;
}

bool ModelBasedIC3::intersects_bad()
{
  solver_->push();

  // assert the last frame (conjunction over clauses)
  assert_frame(reached_k_ + 1);

  // see if it intersects with bad
  solver_->assert_formula(bad_);

  Result r = solver_->check_sat();
  if (r.is_sat()) {
    // create a proof goal for the bad state
    const UnorderedTermSet & statevars = ts_.statevars();
    TermVec cube_vec;
    cube_vec.reserve(statevars.size());
    Term eq;
    for (auto sv : statevars) {
      eq = solver_->make_term(Equal, sv, solver_->get_value(sv));
      cube_vec.push_back(eq);
    }
    Cube c(solver_, cube_vec);
    proof_goals_[reached_k_+1].push_back(c);
  }

  solver_->pop();

  assert(!r.is_unknown());
  return r.is_sat();
}

bool ModelBasedIC3::get_predecessor(size_t i,
                                    const Cube & c,
                                    Cube & out_pred) const
{
  solver_->push();
  assert(i > 0);
  assert(i < frames_.size());
  // F[i-1]
  assert_frame(i - 1);
  // -c
  solver_->assert_formula(solver_->make_term(Not, c.term_));
  // Trans
  solver_->assert_formula(ts_.trans());
  // c'
  solver_->assert_formula(ts_.next(c.term_));

  Result r = solver_->check_sat();
  if (r.is_sat()) {
    out_pred = generalize_predecessor(i, c);
  }

  solver_->pop();
  assert(!r.is_unknown());
  return r.is_sat();
}

ProofGoal ModelBasedIC3::get_next_proof_goal()
{
  assert(has_proof_goals());
  size_t i = proof_goals_.begin()->first;
  ProofGoal pg(proof_goals_.at(i).back(), i);
  // need to remove the proof goal and the map
  // entry if there are no more cubes in the vector
  proof_goals_.at(i).pop_back();
  if (proof_goals_.at(i).empty())
  {
    proof_goals_.erase(i);
  }
  return pg;
}

bool ModelBasedIC3::block_all()
{
  while (has_proof_goals()) {
    ProofGoal pg = get_next_proof_goal();
    // block can fail, which just means a
    // new proof goal will be added
    if (!block(pg) && !pg.second) {
      // if a proof goal cannot be blocked at zero
      // then there's a counterexample
      return false;
    }
  }
  assert(!has_proof_goals());
  return true;
}

bool ModelBasedIC3::block(const ProofGoal & pg)
{
  const Cube & c = pg.first;
  size_t i = pg.second;

  logger.log(
      3, "Attempting to block proof goal <{}, {}>", c.term_->to_string(), i);

  assert(i < frames_.size());
  assert(i >= 0);
  // TODO: assert c -> frames_[i]

  if (i == 0) {
    // can't block anymore -- this is a counterexample
    return false;
  }

  Cube pred;  // populated by get_predecessor if returns false
  if (!get_predecessor(i, c, pred)) {
    // can block this cube
    Term gen_blocking_term = inductive_generalization(i, c);
    logger.log(3, "Blocking term at frame {}: {}", i, c.term_->to_string());
    frames_[i].push_back(gen_blocking_term);
    return true;
  } else {
    // add a new proof goal
    proof_goals_[i - 1].push_back(pred);
    return false;
  }
}

bool ModelBasedIC3::propagate(size_t i)
{
  assert(i + 1 < frames_.size());

  unordered_set<size_t> indices_to_remove;
  TermVec & Fi = frames_.at(i);

  solver_->push();
  assert_frame(i);
  solver_->assert_formula(ts_.trans());

  for (size_t j = 0; j < Fi.size(); ++j) {
    Term &t = Fi[j];

    // Relative inductiveness check
    // Check F[i] /\ t /\ T /\ -t'
    // NOTE: asserting t is redundant because t \in F[i]
    solver_->push();
    solver_->assert_formula(solver_->make_term(Not, ts_.next(t)));

    Result r = solver_->check_sat();
    assert(! r.is_unknown());
    if (r.is_unsat()) {
      // mark for removal
      indices_to_remove.insert(j);
      // add to next frame
      frames_[i + 1].push_back(Fi[j]);
    }

    solver_->pop();
  }

  solver_->pop();

  // keep invariant that terms are kept at highest frame
  // where they are known to hold
  // need to remove some terms from the current frame
  TermVec new_frame_i;
  new_frame_i.reserve(Fi.size() - indices_to_remove.size());
  for (size_t j = 0; j < Fi.size(); ++j) {
    if (indices_to_remove.find(j) == indices_to_remove.end()) {
      new_frame_i.push_back(Fi[j]);
    }
  }

  frames_[i] = new_frame_i;

  return new_frame_i.empty();
}

void ModelBasedIC3::push_frame()
{
  // pushes an empty frame
  frames_.push_back({});
}

Term ModelBasedIC3::inductive_generalization(size_t i, const Cube & c) const
{
  // TODO: implement generalization
  // For now, just getting a blocking clause with no generalization

  // TODO: keep in mind that you cannot drop a literal if it causes c to intersect with the initial states
  assert(!intersects_initial(c.term_));
  return solver_->make_term(Not, c.term_);
}

Clause ModelBasedIC3::down(size_t i, const Clause & c) const
{
  // TODO: implement this when implementing inductive generalization
  // For now, just a stub
  throw PonoException("Not yet implemented");
}

Cube ModelBasedIC3::generalize_predecessor(size_t i, const Cube & c) const
{
  // TODO: do actual generalization
  const UnorderedTermSet & statevars = ts_.statevars();
  TermVec cube_lits;
  cube_lits.reserve(statevars.size());
  for (auto sv : statevars) {
    cube_lits.push_back(solver_->make_term(Equal, sv, solver_->get_value(sv)));
  }
  return Cube(solver_, cube_lits);
}

bool ModelBasedIC3::intersects_initial(const Term & t) const
{
  solver_->push();
  solver_->assert_formula(t);
  solver_->assert_formula(ts_.init());
  Result r = solver_->check_sat();
  solver_->pop();
  return r.is_sat();
}

void ModelBasedIC3::assert_frame(size_t i) const
{
  for (size_t j = i; j < frames_.size(); ++j) {
    for (auto c : frames_[j]) {
      solver_->assert_formula(c);
    }
  }
}

}  // namespace pono
