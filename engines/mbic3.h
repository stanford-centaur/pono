/*********************                                                        */
/*! \file
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
#include <utility>
#include "assert.h"

#include "prover.h"

namespace pono {

// Both clauses and cubes can be represented as vectors of predicates
// They are just negations of eachother

struct Clause
{
  Clause(const smt::SmtSolver & solver, const TermVec & lits) : lits_(lits)
  {
    // sort literals
    std::sort(lits_.begin(), lits_.end(), std::hash<Term>);
    // shouldn't have an empty clause
    assert(lits_.size());

    // create term
    term_ = lits_[0];
    for (size_t i = 1; i < lits_.size(); ++i) {
      term_ = solver->make_term(Or, term_, lits_[i]);
    }
  }
  TermVec lits_;  // list of literals sorted by hash
  Term term_;     // term representation of literals as disjunction
}

struct Cube
{
  Cube(const smt::SmtSolver & solver, const TermVec & lits) : lits_(lits)
  {
    // sort literals
    std::sort(lits_.begin(), lits_.end(), std::hash<Term>);
    // shouldn't have an empty cube
    assert(lits_.size());

    // create term
    term_ = lits_[0];
    for (size_t i = 1; i < lits_.size(); ++i) {
      term_ = solver->make_term(And, term_, lits_[i]);
    }
  }
  TermVec lits_;  // list of literals sorted by hash
  Term term_;     // term representation of literals as conjunction
}

class ModelBasedIC3 : public Prover
{
 public:
  ModelBasedIC3(const Property & p, const smt::SmtSolver & slv);
  ModelBasedIC3(const PonoOptions & opt, const Property p, smt::SmtSolver slv);
  ~ModelBasedIC3();

  typedef Prover super;

  void initialize() override;
  ProverResult check_until(int k) override;

 private:
  /** Check if last frame intersects with bad
   *  @return true iff the last frame intersects with bad
   *  post-condition: if true is returned, bad cube added to proof goals
   */
  bool intersects_bad();
  /** Checks if clause c is relatively inductive to frame i
   *  @param i the frame number
   *  @param c the clause to check
   *  @return true iff c is relatively inductive
   */
  bool rel_ind_check(size_t i, const Clause & c) const;
  /** Gets a new proof goal (and removes it from proof_goals_)
   *  @return a proof goal with the lowest available frame number
   */
  Clause get_next_proof_goal();
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
  bool block(size_t i, const Cube & c);
  /** Propagate all clauses to the highest frame possible */
  void propagate();
  /** Add a new frame*/
  void push_frame();
  /** Attempt to generalize a clause
   *  The standard approach is inductive generalization
   *  @param i the frame number to generalize it against
   *  @param c the clause to generalize
   *  @return a new clause
   */
  Clause generalize_clause(size_t i, const Clause & c) const;
  /** Helper for generalize when using inductive generalization
   *  @param i the frame number
   *  @param c the clause
   *  @return a new clause
   */
  Clause down(size_t i, const Clause & c) const;
  /** Get the predecessor state of a counterexample
   *  @param i the frame number
   *  @param c the bad cube
   *  @return a cube representing the predecessor state
   */
  Cube compute_predecessor(size_t i, const Cube & c) const;
  /** Generalize a counterexample
   *  @param i the frame number
   *  @param c the cube to generalize
   *  @return a new cube
   */
  Cube cex_generalize(size_t i, const Cube & c) const;
  /** Check if the algorithm has found an inductive invariant.
   *  Because clauses are only kept in the highest frame where they still hold
   *  this amounts to checking for any frame that is empty
   *  because this implies that it is equivalent to the next frame
   *  @return true iff the property has been proven
   */
  bool is_proven() const;
  /** Check if the cube intersects with the initial states
   *  @param c the cube to check
   *  @return true iff the cube intersects with the initial states
   */
  bool is_initial(const Cube & c) const;

  // Data structures

  ///< the frames data structure.
  ///< a vector of clauses (except at index 0 which is init)
  ///< for this implementation of IC3
  std::vector<smt::TermVec> frames_;

  ///< outstanding proof goals
  std::unordered_map<size_t, std::vector<Cube>> proof_goals_;

  // useful terms
  Term true_;
  Term false_;
};

}  // namespace pono
