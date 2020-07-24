/*********************                                                  */
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
#include "mbic3.h"

using namespace smt;

namespace pono {

// helpers

static Cube negate(const SmtSolver & slv, const Clause & c)
{
  TermVec cube_lits;
  cube_lits.reserve(c.lits_.size());
  for (auto l : c.lits_) {
    cube_lits.push_back(slv.make_term(Not, l));
  }
  return Cube(slv, cube_lits);
}

static Clause negate(const SmtSolver & slv, const Cube & c)
{
  TermVec clause_lits;
  clause_lits.reserve(c.lits_.size());
  for (auto l : c.lits_) {
    clause_lits.push_back(slv.make_term(Not, l));
  }
  return Clause(slv, clause_lits);
}

ModelBasedIC3::ModelBasedIC3(const Property & p, const smt::SmtSolver & slv)
    : super(p, slv),
      true_(solver_->make_term(true)),
      false_(solver_->make_term(false)),
{
  initialize();
}

ModelBasedIC3::ModelBasedIC3(const PonoOptions & opt,
                             const Property p,
                             smt::SmtSolver slv)
    : super(opt, p, slv),
      true_(solver_->make_term(true)),
      false_(solver_->make_term(false)),
{
  initialize();
}

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
  // TODO: Figure out if we need this
  //       shouldn't have to do this special check
  if (reached_k_ < 1) {
    solver_->push();
    solver_->assert_formula(ts_.init());
    solver_->assert_formula(bad_);
    Result r = solver_->check_sat();
    if (r.is_sat()) {
      return ProverResult::FALSE;
    } else {
      assert(r.is_unsat());
      reached_k_ = 1;  // keep reached_k_ aligned with number of frames
    }
    solver_->pop();
  }

  while (reached_k_ <= k) {
    assert(reached_k_ == frames_.size());
    // blocking phase
    while (intersects_bad()) {
      assert(!proof_goals_.empty());
      if (!block_all()) {
        // counter-example
        return ProverResult::FALSE;
      }
    }

    // propagation phase
    push_frame();
    for (size_t i = 1; i < frames_.size() - 1; ++i) {
      if (propagate(i)) {
        return ProverResult::TRUE;
      }
    }
  }

  return ProverResult::UNKNOWN;
}

bool ModelBasedIC3::intersects_bad()
{
  solver_->push();

  // assert the frame (conjunction over clauses)
  assert_frame(reached_k_);

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
    proof_goals_[reached_k_].push_back(c);
  }

  solver_->pop();

  assert(!r.is_unknown());
  return r.is_sat();
}

bool ModelBasedIC3::get_predecessor(size_t i, const Cube & c, Cube & out_cti)
{
  solver_->push();
  assert(i < frames_.size());
  // F[i]
  assert_frame(i);
  // -c
  solver_->assert_formula(solver_->make_term(Not, c.term_));
  // Trans
  solver_->assert_formula(ts_.trans());
  // c'
  solver_->assert_formula(ts_.next(c.term_));

  Result r = solver_->check_sat();
  if (r.is_sat()) {
    const UnorderedTermSet & statevars = ts_.statevars();
    TermVec cube_lits;
    cube_lits.reserve(statevars.size());
    for (auto sv : statevars) {
      cube_lits.push_back(
          solver_->make_term(Equal, sv, solver_->get_value(sv)));
    }
    out_cti = Cube(solver_, cube_lits);
  }

  solver_->pop();
  assert(!r.is_unknown());
  return r.is_unsat();
}

ProofGoal ModelBasedIC3::get_next_proof_goal()
{
  assert(has_proof_goals());
  size_t i = proof_goals_.begin()->first;
  ProofGoal pg(proof_goals_.at(i).back(), i);
  // need to remove the proof goal
  proof_goals_.at(i).pop_back();
  return pg;
}

bool ModelBasedIC3::block_all()
{
  while (has_proof_goals()) {
    ProofGoal pg = get_next_proof_goal();
    // block can fail, which just means a
    // new proof goal will be added
    if (pg.first && !block(pg)) {
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
  size_t i = pg.first;
  Cube & c = pg.second;

  assert(i < frames_.size());
  assert(i >= 0);

  if (i == 0) {
    // can't block anymore -- this is a counterexample
    return false;
  }

  Cube cti;  // populated by get_predecessor if returns false
  if (get_predecessor(i-1, c, cti)) {
    // can block this cube
    Clause neg_c = negate(solver_, c);
    Clause gen_blocking_clause = generalize_clause(i, neg_c);
    frames_[i].push_back(gen_blocking_clause);
  } else {
    // add a new proof goal
    proof_goals_[i - 1].push_back(generalize_cti(i - 1, cti));
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
    solver_->assert_formula(solver_->make_term(Not, ts_next(t)));

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

Clause ModelBasedIC3::generalize_clause(size_t i, const Clause & c) const
{
  // TODO: actual generalization
  // For now, just a NOP stub
  return c;
}

Clause ModelBasedIC3::down(size_t i, const Clause & c) const
{
  // TODO: implement this when implementing inductive generalization
  // For now, just a stub
  throw CosaException("Not yet implemented");
}

Cube ModelBasedIC3::generalize_cti(size_t i, const Cube & c) const
{
  // TODO: implement this
  // For now, just a NOP stub.
  return c;
}

bool ModelBasedIC3::is_initial(const Cube & c) const
{
  solver_->push();
  solver_->assert_formula(c.term_);
  solver_->assert_formula(ts_.init());
  Result r = solver_->check_sat();
  solver_->pop();
  return r.is_sat();
}

void ModelBasedIC3::assert_frame(size_t i) const
{
  for (size_t j = i; j <= reached_k_; ++j) {
    for (auto c : frames_[j]) {
      solver_->assert_formula(c);
    }
  }
}

}  // namespace pono
