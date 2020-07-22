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

using ProofGoal = std::pair<Cube, size_t>;

class ModelBasedIC3 : public Prover
{
 public:
  IC3(const Property & p, smt::SmtSolver slv);
  IC3(const PonoOptions & opt, const Property p, smt::SmtSolver slv);
  ~IC3();
  ProverResult prove() override;
  ProverResult check_until(int k) override;
  void initialize() override;

 private:
  /** Checks if clause c is relatively inductive to frame i
   *  @param i the frame number
   *  @param c the clause to check
   *  @return true iff c is relatively inductive
   */
  bool rel_ind_check(size_t i, Clause c);
  /** Gets a new proof goal
   *  @return a proof goal with the lowest available frame number
   */
  Clause get_next_proof_goal();
  /** Attempt to block cube c at frame i
   *  @param i the frame number
   *  @param c the cube to try blocking
   *  @return true iff the cube was blocked, otherwise a new proof goal was
   * added to the proof goals
   */
  bool block(size_t i, Cube c);
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
  Clause generalize_clause(size_t i, Clause c);
  /** Helper for generalize when using inductive generalization
   *  @param i the frame number
   *  @param c the clause
   *  @return a new clause
   */
  Clause down(size_t i, Clause c);
  /** Get the predecessor state of a counterexample
   *  @param i the frame number
   *  @param c the bad cube
   *  @return a cube representing the predecessor state
   */
  Cube compute_predecessor(size_t i, Cube c);
  /** Generalize a counterexample
   *  @param i the frame number
   *  @param c the cube to generalize
   *  @return a new cube
   */
  Cube cex_generalize(size_t i, Cube c);
  /** Check if the algorithm has found an inductive invariant.
   *  Because clauses are only kept in the highest frame where they still hold
   *  this amounts to checking for any frame that is empty
   *  because this implies that it is equivalent to the next frame
   *  @return true iff the property has been proven
   */
  bool is_proven();

  // Data structures

  ///< the frames data structure
  std::vector<std::vector<Clause>> frames_;

  ///< outstanding proof goals
  std::unordered_map<size_t, std::vector<Cube>> proof_goals_;
};

}  // namespace pono
