/*********************                                                  */
/*! \file mbic3.h
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
#pragma once

#include <algorithm>
#include <map>
#include <utility>
#include "assert.h"

#include "prover.h"

namespace pono {

// Both clauses and cubes can be represented as vectors of predicates
// They are just negations of eachother

/** Less than comparison the hash of two terms
 *  for use in sorting
 *  @param t0 the first term
 *  @param t1 the second term
 *  @return true iff t0's hash is less than t1's hash
 */
bool term_hash_lt(const smt::Term & t0, const smt::Term & t1)
{
  return (t0->hash() < t1->hash());
}

struct Clause
{
  Clause() {}
  Clause(const smt::SmtSolver & solver, const smt::TermVec & lits) : lits_(lits)
  {
    // sort literals
    std::sort(lits_.begin(), lits_.end(), term_hash_lt);
    // shouldn't have an empty clause
    assert(lits_.size());

    // create term
    term_ = lits_[0];
    for (size_t i = 1; i < lits_.size(); ++i) {
      term_ = solver->make_term(smt::Or, term_, lits_[i]);
    }
  }

  smt::TermVec lits_;  // list of literals sorted by hash
  smt::Term term_;     // term representation of literals as disjunction
};

struct Cube
{
  Cube() {}
  Cube(const smt::SmtSolver & solver, const smt::TermVec & lits) : lits_(lits)
  {
    // sort literals
    std::sort(lits_.begin(), lits_.end(), term_hash_lt);
    // shouldn't have an empty cube
    assert(lits_.size());

    // create term
    term_ = lits_[0];
    for (size_t i = 1; i < lits_.size(); ++i) {
      term_ = solver->make_term(smt::And, term_, lits_[i]);
    }
  }
  smt::TermVec lits_;  // list of literals sorted by hash
  smt::Term term_;     // term representation of literals as conjunction
};

using ProofGoal = std::pair<Cube, size_t>;

class ModelBasedIC3 : public Prover
{
 public:
  // TODO: make references const once engine-fresh-solver branch is merged
  ModelBasedIC3(const Property & p, smt::SmtSolver & slv);
  ModelBasedIC3(const PonoOptions & opt, const Property p, smt::SmtSolver slv);
  virtual ~ModelBasedIC3();

  typedef Prover super;

  void initialize() override;
  ProverResult check_until(int k) override;

 private:
  /** Perform a IC3 step
   *  @param i
   */
  ProverResult step(int i);
  /** Perform the base IC3 step (zero case)
   */
  ProverResult step_0();

  /** Check if last frame intersects with bad
   *  @return true iff the last frame intersects with bad
   *  post-condition: if true is returned, bad cube added to proof goals
   */
  bool intersects_bad();
  /** Get the predecessor of a cube c in frame i
   *  aka see if c is reachable from frame i-1
   *  @requires c -> F[i]
   *  @param i the frame number
   *  @param t the term to check
   *  @param cti the cube to populate with a cti
   *  @return true iff c is reachable from the frame i
   *  @ensures returns true  : pred -> F[i-1] /\ (pred, c) \in [T]
   *           returns false : pred unchanged, F[i-1] /\ T /\ c' is unsat
   */
  bool get_predecessor(size_t i, const Cube & c, Cube & out_pred) const;
  /** Check if there are more proof goals
   *  @return true iff there are more proof goals
   */
  bool has_proof_goals() const { return !proof_goals_.empty(); };
  /** Gets a new proof goal (and removes it from proof_goals_)
   *  @requires has_proof_goals()
   *  @return a proof goal with the lowest available frame number
   *  @alters proof_goals_
   *  @ensures returned proof goal is from lowest frame in proof goals
   */
  ProofGoal get_next_proof_goal();
  /** Attempt to block all proof goals
   *  to ensure termination, always choose proof goal with
   *  smallest time
   *  @return true iff all proof goals were blocked
   */
  bool block_all();
  /** Attempt to block cube c at frame i
   *  @param i the frame number
   *  @param c the cube to try blocking
   *  @return true iff the cube was blocked, otherwise a new proof goal was
   * added to the proof goals
   */
  bool block(const ProofGoal & pg);
  /** Try propagating all clauses from frame index i to the next frame.
   *  @param i the frame index to propagate
   *  @return true iff all the clauses are propagated (this means property was
   * proven)
   */
  bool propagate(size_t i);
  /** Add a new frame*/
  void push_frame();
  /** Attempt to generalize a clause
   *  The standard approach is inductive generalization
   *  @requires !get_predecessor(i, c, _)
   *  @param i the frame number to generalize it against
   *  @param c the cube to find a general predecessor for
   *  @return a new term P
   *  @ensures P -> !c /\ F[i-1] /\ P /\ T /\ !P' is unsat
   */
  smt::Term inductive_generalization(size_t i, const Cube & c) const;
  /** Helper for generalize when using inductive generalization
   *  @param i the frame number
   *  @param c the clause
   *  @return a new clause
   */
  Clause down(size_t i, const Clause & c) const;
  /** Generalize a counterexample
   *  @requires get_predecessor(i, c)
   *  @param i the frame number
   *  @param c the cube to generalize
   *  @return a new cube d
   *  @ensures d -> F[i-1] /\ forall s \in [d] exists s' \in [c]. (d,c) \in [T]
   */
  Cube generalize_predecessor(size_t i, const Cube & c) const;
  /** Check if the term intersects with the initial states
   *  @param t the term to check
   *  @return true iff t intersects with the initial states
   */
  bool intersects_initial(const smt::Term & t) const;
  /** Add all the terms at Frame i
   *  Note: the frames_ data structure keeps terms only in the
   *  highest frame where they are known to hold
   *  Thus, asserting Fi actually needs to add terms
   *  from Fi and all frames after it
   *  @param i the frame number
   */
  void assert_frame(size_t i) const;

  // Data structures

  ///< the frames data structure.
  ///< a vector of clauses (except at index 0 which is init)
  ///< for this implementation of IC3
  std::vector<smt::TermVec> frames_;

  ///< outstanding proof goals -- kept sorted so lower ones
  ///< can be handled first
  std::map<size_t, std::vector<Cube>> proof_goals_;

  // useful terms
  smt::Term true_;
  smt::Term false_;
};

}  // namespace pono
