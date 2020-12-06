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
** \brief Simple implementation of IC3 using model values
**
**/
#pragma once

#include <algorithm>
#include <map>
#include <utility>
#include "assert.h"

#include "smt-switch/term_translator.h"

#include "engines/prover.h"

namespace pono {

/** Represents a conjunction of boolean terms
 *  Keeps the conjuncts in a vector as well
 *  and sorts them
 */
struct Conjunction
{
  Conjunction() {}
  Conjunction(const smt::SmtSolver & solver, const smt::TermVec & conjuncts);
  smt::TermVec conjuncts_;  // list of literals sorted by hash
  smt::Term term_;          // term representation of literals as conjunction
};

struct ProofGoal
{
  // based on open-source ic3ia ProofObligation
  Conjunction conj;
  size_t idx;
  // TODO: refactor this. shared_ptr probably overkill
  std::shared_ptr<ProofGoal> next;

  ProofGoal(Conjunction c, size_t i, std::shared_ptr<ProofGoal> n)
      : conj(c), idx(i), next(n)
  {
  }
};

class ModelBasedIC3 : public Prover
{
 public:
  ModelBasedIC3(Property & p, smt::SolverEnum se);
  ModelBasedIC3(Property & p, const smt::SmtSolver & slv);
  ModelBasedIC3(const PonoOptions & opt,
                Property & p,
                smt::SolverEnum se);
  ModelBasedIC3(const PonoOptions & opt,
                Property & p,
                const smt::SmtSolver & slv);
  virtual ~ModelBasedIC3();

  typedef Prover super;

  void initialize() override;

  /** Set up the interpolator_
   */
  void initialize_interpolator();

  ProverResult check_until(int k) override;
  bool witness(std::vector<smt::UnorderedTermMap> & out) override;

 protected:
  /** Perform a IC3 step
   *  @param i
   */
  ProverResult step(int i);
  /** Perform the base IC3 step (zero case)
   */
  ProverResult step_0();

  /** Check if a transition from the second to last frame can result in a bad
   * state
   *  @return true iff the last frame intersects with bad
   *  post-condition: if true is returned, bad cube added to proof goals
   */
  virtual bool intersects_bad();
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
  virtual bool get_predecessor(size_t i,
                               const Conjunction & c,
                               Conjunction & out_pred);
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
  /** Create and add a proof goal for cube c for frame i
   *  @param c the cube of the proof goal
   *  @param i the frame number for the proof goal
   *  @param n pointer to the proof goal that led to this one -- null for bad
   *  (i.e. end of trace)
   */
  void add_proof_goal(const Conjunction & c,
                      size_t i,
                      std::shared_ptr<ProofGoal> n);
  /** Attempt to block all proof goals
   *  to ensure termination, always choose proof goal with
   *  smallest time
   *  @return true iff all proof goals were blocked
   */
  virtual bool block_all();
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
  /** Add a new frame */
  void push_frame();
  /** Adds a constraint to frame i and (implicitly) all frames below it
   *  @param i highest frame to add constraint to
   *  @param constraint the constraint to add
   */
  void constrain_frame(size_t i, const smt::Term & constraint);
  /** Attempt to generalize a clause
   *  The standard approach is inductive generalization
   *  @requires !get_predecessor(i, c, _)
   *  @param i the frame number to generalize it against
   *  @param c the cube to find a general predecessor for
   *  @return a new term P
   *  @ensures P -> !c /\ F[i-1] /\ P /\ T /\ !P' is unsat
   */
  smt::Term inductive_generalization(size_t i, const Conjunction & c);
  /** Generalize a counterexample
   *  @requires get_predecessor(i, c)
   *  @param i the frame number
   *  @param c the cube to generalize
   *  @return a new cube d
   *  @ensures d -> F[i-1] /\ forall s \in [d] exists s' \in [c]. (d,c) \in [T]
   */
  Conjunction generalize_predecessor(size_t i, const Conjunction & c);
  /** Check if there are common assignments
   *  between A and B
   *  i.e. if A /\ B is SAT
   *  @param A the first term
   *  @param B the second term
   *  @return true iff there is an intersection
   */
  bool intersects(const smt::Term & A, const smt::Term & B);
  /** Check if the term intersects with the initial states
   *  syntactic sugar for intersects(ts_->init(), t);
   *  @param t the term to check
   *  @return true iff t intersects with the initial states
   */
  bool intersects_initial(const smt::Term & t);
  /** Add all the terms at Frame i
   *  Note: the frames_ data structure keeps terms only in the
   *  highest frame where they are known to hold
   *  Thus, asserting Fi actually needs to add terms
   *  from Fi and all frames after it
   *  @param i the frame number
   */
  void assert_frame_labels(size_t i) const;

  smt::Term get_frame(size_t i) const;

  void assert_trans_label() const;

  virtual smt::Term get_trans() const;

  void fix_if_intersects_initial(smt::TermVec & to_keep,
                                 const smt::TermVec & rem);

  size_t push_blocking_clause(size_t i, smt::Term c);

  /** Helper function to reduce assumptions using unsat cores.
   *  @param input formula
   *  @param vector of assumptions
   *  @param vector to store reduced assumptions
   *  @param vector to store removed assumptions (if not NULL)
   */
  void reduce_assump_unsatcore(const smt::Term &formula,
                               const smt::TermVec &assump,
                               smt::TermVec &out_red,
                               smt::TermVec *out_rem = NULL);

  smt::Term label(const smt::Term & t);

  /** Creates a reduce and of the vector of boolean terms
   *  It also sorts the vector by the hash
   *  Note: this will fail for empty vectors
   *  @param vec the vector of boolean terms
   *  @param slv (optional) the solver to use, defaults to solver_
   *  @return the conjunction of all the terms
   */
  smt::Term make_and(smt::TermVec vec, smt::SmtSolver slv = NULL) const;

  /** Creates a reduce or of the vector of boolean terms
   *  It also sorts the vector by the hash
   *  Note: this will fail for empty vectors
   *  @param vec the vector of boolean terms
   *  @param slv (optional) the solver to use, defaults to solver_
   *  @return the disjunction of all the terms
   */
  smt::Term make_or(smt::TermVec vec, smt::SmtSolver slv = NULL) const;

  /** Pushes a solver context and keeps track of the context level
   *  updates solver_context_
   */
  void push_solver_context();

  /** Pops a solver context and keeps track of the context level
   *  updates solver_context_
   */
  void pop_solver_context();

  /** Sets the semantics for labels (if not already set)
   *  e.g. init_label_, and trans_label_
   *  -- checks if they are null and if so, gives semantics
   *  Note: this might be overloaded in derived classes
   *  Thus, for now should not be called in the constructor
   *  TODO see if there is a cleaner way to do this
   *       want to share the labels between this class and derived
   *       but change the associated transition relation
   *       need to be careful with calls to virtual functions
   *       in the constructor
   */
  virtual void set_labels();

  /** Analyze TS to see if there are any unsupported logics
   */
  virtual void check_ts() const;

  /** Check that term is only over current state variables
   *  This method is virtual so that derived classes can
   *  use a different transition system (e.g. an abstract one)
   *  for the check
   *  @param t the term to check
   *  @return true iff the term only containts current state variables
   */
  virtual bool only_curr(smt::Term & t);

  /** Map current state to next state variables
   *  This automatically uses the correct transition system
   *  i.e. in derived classes operating over an abstract system, it
   *  will use the abstraction
   */
  virtual smt::Term next(const smt::Term & t) const;

  /** Sets the invariant with the contents of frame i
   *  @param i the frame holding the invariant
   *  Note: some IC3 implementations might need to do some translation
   *  from an abstraction
   */
  void set_invar(size_t i);

  // Data structures

  ///< the frames data structure.
  ///< a vector of terms (typically clauses but depends on mode)
  ///< for this implementation of IC3
  std::vector<smt::TermVec> frames_;

  // labels for activating assertions
  smt::Term init_label_;          ///< label to activate init
  smt::Term trans_label_;         ///< label to activate trans
  smt::TermVec frame_labels_;     ///< labels to activate frames
  smt::UnorderedTermMap labels_;  //< labels for unsat cores

  ///< stack of outstanding proof goals
  std::vector<ProofGoal> proof_goals_;

  // useful terms from the default solver
  smt::Term solver_true_;
  smt::Term solver_false_;

  // NOTE: need to be sure to always use [push|pop]_solver_context
  // instead of doing it directly on the solver or this will be
  // out of sync
  // TODO could consider adding this to smt-switch
  //      but would have to be implemented in each solver
  size_t solver_context_;  ///< the current context of the solver

  // for ic3_indgen_mode_ == 2
  // interpolant based generalization
  smt::SmtSolver interpolator_;
  std::unique_ptr<smt::TermTranslator> to_interpolator_;
  std::unique_ptr<smt::TermTranslator> to_solver_;
};

class DisjointSet
{
 public:
  DisjointSet();
  ~DisjointSet();

  void add(const smt::Term & a, const smt::Term & b);
  smt::Term find(const smt::Term & t);
  void clear();

 private:
  // member to group's leader
  smt::UnorderedTermMap leader_;
  // group leader to group
  std::unordered_map<smt::Term, smt::UnorderedTermSet> group_;
};

}  // namespace pono
