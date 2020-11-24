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

#include <type_traits>

#include "assert.h"
#include "smt/available_solvers.h"

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

/** IC3Base */

IC3Base::IC3Base(Property & p, smt::SolverEnum se, IC3UnitCreator ic)
    : super(p, se), mk_unit(ic), reducer_(create_solver(se)), solver_context_(0)
{
  initialize();
}

IC3Base::IC3Base(Property & p, const smt::SmtSolver & s, IC3UnitCreator ic)
    : super(p, s),
      mk_unit(ic),
      reducer_(create_solver(s->get_solver_enum())),
      solver_context_(0)
{
  initialize();
}

IC3Base::IC3Base(const PonoOptions & opt,
                 Property & p,
                 smt::SolverEnum se,
                 IC3UnitCreator ic)
    : super(opt, p, se),
      mk_unit(ic),
      reducer_(create_solver(se)),
      solver_context_(0)
{
  initialize();
}
IC3Base::IC3Base(const PonoOptions & opt,
                 Property & p,
                 const smt::SmtSolver & s,
                 IC3UnitCreator ic)
    : super(opt, p, s),
      mk_unit(ic),
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

  solver_true_ = solver_->make_term(true);
}

ProverResult IC3Base::check_until(int k) { throw PonoException("NYI"); }

bool IC3Base::witness(std::vector<smt::UnorderedTermMap> & out)
{
  throw PonoException("NYI");
}

// Protected Methods

bool IC3Base::intersects_bad() { throw PonoException("NYI"); }

ProverResult IC3Base::step(int i) { throw PonoException("NYI"); }

ProverResult IC3Base::step_0() { throw PonoException("NYI"); }

// Helper methods

bool IC3Base::block_all() { throw PonoException("NYI"); }

bool IC3Base::block(const IC3Goal & pg) { throw PonoException("NYI"); }

bool IC3Base::propagate(size_t i)
{
  assert(i + 1 < frames_.size());

  unordered_set<size_t> indices_to_remove;
  const vector<IC3Unit> & Fi = frames_.at(i);

  push_solver_context();
  assert_frame_labels(i);
  assert_trans_label();

  for (size_t j = 0; j < Fi.size(); ++j) {
    const Term & t = Fi.at(j).get_term();

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
  assert(i < frame_labels_.size());
  assert(frame_labels_.size() == frames_.size());
  solver_->assert_formula(
      solver_->make_term(Implies, frame_labels_.at(i), constraint.get_term()));
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
      res = solver_->make_term(And, res, u.get_term());
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
  IC3Goal pg = std::move(proof_goals_.back());
  proof_goals_.pop_back();
  return pg;
}

void IC3Base::add_proof_goal(const IC3Unit & c, size_t i, unique_ptr<IC3Goal> n)
{
  // IC3Unit aligned with frame so proof goal should be negated
  // e.g. for bit-level IC3, IC3Unit is a Clause and the proof
  // goal should be a Cube
  assert(c.is_negated());
  proof_goals_.push_back(IC3Goal(c, i, std::move(n)));
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
  assert(!u.is_negated());
  Term c = u.get_term();
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
