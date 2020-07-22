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

Cube negate(const SmtSolver & slv, const Clause & c)
{
  TermVec cube_lits;
  cube_lits.reserve(c.lits_.size());
  for (auto l : c.lits_) {
    cube_lits.push_back(slv.make_term(Not, l));
  }
  return Cube(slv, cube_lits);
}

Clause negate(const SmtSolver & slv, const Cube & c)
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
    propagate();
    if (is_proven()) {
      return ProverResult::TRUE;
    }
  }

  return ProverResult::UNKNOWN;
}

bool ModelBasedIC3::intersects_bad()
{
  solver_->push();

  // assert the frame (conjunction over clauses)
  for (auto c : frames_.back()) {
    solver_->assert_formula(c);
  }

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

bool ModelBasedIC3::rel_ind_check(size_t i, const Clause & c)
{
  solver_->push();
  Result r = rel_ind_check_helper(i, c);
  solver_->pop();

  assert(!r.is_unknown());
  return r.is_unsat();
}

bool ModelBasedIC3::rel_ind_check(size_t i, const Clause & c, Cube & cti)
{
  solver_->push();
  Result r = rel_ind_check_helper(i, c);

  if (r.is_sat()) {
    const UnorderedTermSet & statevars = ts_.statevars();
    TermVec cube_lits;
    cube_lits.reserve(statevars.size());
    for (auto sv : statevars) {
      cube_lits.push_back(
          solver_->make_term(Equal, sv, solver_->get_value(sv)));
    }
    cti = Cube(solver_, cube_lits);
  }

  solver_->pop();
  assert(!r.is_unknown());
  return r.is_unsat();
}

Result ModelBasedIC3::rel_ind_check_helper(size_t i, const Clause & c)
{
  // Check F[i] /\ -c /\ T /\ c'
  assert(i < frames_.size());

  // F[i]
  for (auto c : frames_[i]) {
    solver_->assert_formula(c);
  }

  // -c
  Cube neg_c = negate(solver_, c);
  solver_->assert_formula(neg_c.term_);

  // T
  solver_->assert_formula(ts_.trans());

  // c'
  Term cprime = ts_.next(c.term_);
  solver_->assert_formula(cprime);

  return solver_->check_sat();
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

bool block_all()
{
  while (has_proof_goals()) {
    ProofGoal pg = get_next_proof_goal();
    // block can fail, which just means a
    // new proof goal will be added
    if (!block(pg) && pg.first == 0) {
      // if a proof goal cannot be blocked at zero
      // then there's a counterexample
      return false;
    }
  }
  assert(!has_proof_goals());
  return true;
}

}  // namespace pono
