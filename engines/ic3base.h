/*********************                                                  */
/*! \file ic3base.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan
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
**        To create a particular IC3 instantiation, you must implement the
*following:
**           - implement get_model_ic3_formula and give it semantics to produce
**             the corresponding IC3Formula for your flavor of IC3
**             (assumes solver_'s state is SAT from a failed rel_ind_check)
**           - optionally override inductive_generalization
**             (default aggressively minimizes)
**           - optionally override generalize_predecessor
**             (default just takes whole model)
**             this one is special because it is called with solver_context == 1
**             and assuming the state is SAT so you can query with
**             solver_->get_value(Term)
**             However, this also means you can't use the solver_ because
**             it's polluted with other assertions. All calls are expected
**             to be through a separate solver (e.g. through the reducer_)
**           - implement check_ts which just checks if there are any
**             theories / syntax used in the transition system which
**             is not supported by this instantiation
**          [OPTIONAL]
**           - implement all the ic3formula_* functions for creating and
**             manipulating an IC3Formula if the defaults are not right
**           - implement abstract() and refine() if this is a CEGAR
**             flavor of IC3
**           - override reset_solver if you need to add back in constraints
**             to the reset solver that aren't handled by the default
**             implementation
**           - if you decide to use additional labels, make sure to add
**             them to solver_ in both initialize and reset_solver
**           - if your implementation tends to use the same labels often
**             you can add the semantics (e.g. equality / implication)
**             at the base context level, and just make sure that
**             is_global_label returns true for those labels by overriding
**             just don't forget to re-add those assertions in reset_solver
**           - if your algorithm has approximate predecessor generalization
**             meaning it is not guaranteed to keep the predecessor in F[i-1]
**             without intersecting with a lower frame (e.g. F[i-2]), then
**             set approx_pregen_ to true and the base class will handle
**             maintaining that invariant of IC3
**
**        Important Notes:
**           - be sure to use [push/pop]_solver_context instead of using
**             the solver interface directly so that the assertions on
**             the solver context (tracked externally) are correct
**           - be sure to use the check_sat() and check_sat_assuming
**             member methods -- they will keep track of the number of
**             calls since a reset
**
**/
#pragma once

#include <algorithm>
#include <cstddef>
#include <queue>
#include <vector>

#include "core/proverresult.h"
#include "core/refineresult.h"
#include "core/ts.h"
#include "engines/prover.h"
#include "options/options.h"
#include "smt-switch/smt.h"
#include "smt-switch/utils.h"

namespace pono {

struct IC3Formula
{
  // nullary constructor
  IC3Formula() : term(nullptr) {}
  IC3Formula(const smt::Term & t, const smt::TermVec & c, bool n)
      : term(t), children(c), disjunction(n)
  {
    std::sort(children.begin(), children.end());
  }

  IC3Formula(const IC3Formula & other)
      : term(other.term),
        children(other.children),
        disjunction(other.disjunction)
  {
  }

  virtual ~IC3Formula() {}

  /** Returns true iff this IC3Formula has not been initialized */
  bool is_null() const { return (term == nullptr); };

  smt::Term term;  ///< term representation of this formula
  smt::TermVec
      children;  ///< flattened children of either a disjunction or conjunction
  bool disjunction;  ///< true if currently representing a disjunction
  // NOTE: treating the disjunction as the primary orientation
  //       e.g. aligned with what's kept in frames (disjunctions), not proof
  //       goals (conjunctions), but IC3Formula can represent both
};

struct ProofGoal
{
  // based on open-source ic3ia ProofObligation
  IC3Formula target;
  size_t idx;
  const ProofGoal * next;

  ProofGoal(IC3Formula u, size_t i, const ProofGoal * n)
      : target(u), idx(i), next(n)
  {
  }
};

/**
 * Ordering for proof obligations in the priority queue (see below) -- borrowed
 * from open-source ic3ia implementation
 */
struct ProofGoalOrder
{
  // comparison for priority queue
  // since priority queue returns largest element, we swap the arguments
  // -- we want the lowest index to be processed first
  bool operator()(const ProofGoal * a, const ProofGoal * b) const
  {
    return b->idx < a->idx;
  }
};

/**
 * Priority queue of proof obligations inspired by open-source ic3ia
 * implementation
 */
class ProofGoalQueue
{
 public:
  ~ProofGoalQueue();

  void clear();
  void new_proof_goal(const IC3Formula & c,
                      unsigned int t,
                      const ProofGoal * n = NULL);
  ProofGoal * top();
  void pop();
  bool empty() const;

 private:
  std::priority_queue<ProofGoal *, std::vector<ProofGoal *>, ProofGoalOrder>
      queue_;
  std::vector<ProofGoal *> store_;
};

class IC3Base : public SafetyProver
{
 public:
  typedef SafetyProver super;

  /** IC3Base constructors take the normal arguments for a prover
   *  + a function that can create an IC3Formula
   *  Depending on the derived class IC3 implementation, the exact
   *  type of IC3Formula will differ: e.g. Clause, Disjunction
   */
  IC3Base(const SafetyProperty & p,
          const TransitionSystem & ts,
          const smt::SmtSolver & s,
          PonoOptions opt = PonoOptions());

  virtual ~IC3Base();

  void initialize() override;

  ProverResult check_until(int k) override;

  size_t witness_length() const override;

 protected:
  bool compute_witness() override;

  bool compute_witness(const TransitionSystem & ts);

  smt::UnsatCoreReducer reducer_;

  ///< keeps track of the current context-level of the solver
  // NOTE: if solver is passed in, it could be off
  //       currently no way to check
  //       then solver_context_ is relative to the starting context
  size_t solver_context_;

  size_t num_check_sat_since_reset_;

  bool failed_to_reset_solver_;  ///< some solvers don't support reset
                                 ///< assertions. Stop trying for those solvers.

  smt::TermVec cex_;  ///< a vector of terms over state variables describing
                      ///< a (possibly abstract) counterexample trace

  bool approx_pregen_;  ///< if set to true then predecessor generalization
                        ///< might over-generalize leading to intersection
                        ///< with F[i-2]
                        ///< in this case, it will do an extra call to fix
                        ///< this to maintain the invariant that the predecessor
                        ///< is within F[i-1] but not F[i-2]
                        ///< default is false and should be set in algorithms
                        ///< with approximate predecessor generalization

  ///< the frames data structure.
  ///< a vector of the given Unit template
  ///< which changes depending on the implementation
  std::vector<std::vector<IC3Formula>> frames_;

  ///< priority queue of outstanding proof goals
  // labels for activating assertions
  smt::Term init_label_;          ///< label to activate init
  smt::Term trans_label_;         ///< label to activate trans
  smt::Term bad_label_;           ///< label to activate bad
  smt::TermVec frame_labels_;     ///< labels to activate frames
  smt::UnorderedTermMap labels_;  ///< labels for unsat cores

  // useful terms
  smt::Term solver_true_;

  smt::Sort boolsort_;

  // re-usable data structures
  // NOTE: be sure not to overwrite these in nested function calls
  smt::TermVec assumps_;  ///< used for storing assumptions

  // TODO Make sure all comments are updated!

  // *************************** Main Methods *********************************

  // ************************** Virtual Methods *******************************
  // IMPORTANT for derived classes
  // These methods should be implemented by a derived class for a particular
  // "flavor" of IC3

  // Pure virtual methods that must be overridden

  /** Get an IC3Formula from the current model
   *  @requires last call to check_sat of solver_ was satisfiable and context
   * hasn't changed
   *  @return an IC3Formula over current state variables with is_disjunction
   *          false depending on the flavor of IC3, this might be a boolean
   *          cube, a theory cube, a cube of predicates, etc...
   */
  virtual IC3Formula get_model_ic3formula() const = 0;

  /** Check whether a given IC3Formula is valid
   *  e.g. if this is a boolean clause it would
   *    check that it's a disjunction of literals
   *  (for debugging)
   *  @param u the IC3Formula to check
   *  @return true iff this is a valid IC3Formula for this
   *          flavor of IC3
   */
  virtual bool ic3formula_check_valid(const IC3Formula & u) const = 0;

  /** Checks if every thing in the current transition system is supported
   *  by the current instantiation
   *  throws a PonoException with a relevant message if not.
   */
  virtual void check_ts() const = 0;

  // virtual methods that can optionally be overridden

  /** Creates a disjunction IC3Formula from a vector of terms
   *  @param c the children terms
   *  @ensures resulting IC3Formula children == c
   *  @ensures resulting IC3Formula with is_disjunction true
   */
  virtual IC3Formula ic3formula_disjunction(const smt::TermVec & c) const;

  /** Creates a conjunction IC3Formula from a vector of terms
   *  @param c the children terms
   *  @ensures resulting IC3Formula children == c
   *  @ensures resulting IC3Formula with is_disjunction false
   *  e.g. for a ClauseHandler, this method will create a cube
   *  note: assumes the children are already in the right polarity
   *  (doesn't negate them)
   */
  virtual IC3Formula ic3formula_conjunction(const smt::TermVec & c) const;

  /** Negates an IC3Formula
   *  @param u the IC3Formula to negate
   */
  virtual IC3Formula ic3formula_negate(const IC3Formula & u) const;

  /** Attempt to generalize before blocking in frame i
   *  The standard approach is inductive generalization
   *  @requires !rel_ind_check(i, c, _)
   *  @requires c is a conjunction (e.g. !c.disjunction)
   *  @param i the frame number to generalize it against
   *  @param c the IC3Formula that should be blocked
   *  @return a generalized IC3Formula
   *
   *  Let the returned formula be d
   *  @ensures d -> !c and F[i-1] /\ d /\ T /\ !d' is unsat
   *           e.g. it blocks c and is inductive relative to F[i-1]
   */
  virtual IC3Formula inductive_generalization(size_t i, const IC3Formula & c);

  /** Generalize a counterexample
   *  @requires rel_ind_check(i, c)
   *  @requires the solver_ context is currently satisfiable
   *  @param i the frame number
   *  @param c the term that failed to be blocked at i
   *         (over current state vars)
   *  @param pred the predecessor. originally passed as full assignment
   *         and then is updated in to be a generalized predecessor
   *  @ensures pred -> F[i-1] /\
   *           if !approx_pregen_ then (F[i-2] /\ pred) is unsat
   *           Most algorithms also ensure:
   *           forall s \in [pred] exists s' \in [c]. (pred ,c) \in [T]
   *  @ensures no calls to the solver_ because the context is polluted with
   *           other assertions
   */
  virtual void predecessor_generalization(size_t i,
                                          const smt::Term & c,
                                          IC3Formula & pred);

  /** Generates an abstract transition system
   *  Typically this would set the ts_ pointer to the abstraction
   *  so that the algorithm can run on the abstraction
   *
   *  by default this does nothing, e.g. meant for an algorithm that does not
   *  do abstraction refinement
   */
  virtual void abstract()
  {
    // by default this is a No-Op
  }

  /** Refines an abstract transition system
   *  By default, this does nothing (e.g. assumes ts_ is not abstract)
   *  Any CEGAR implementations of IC3 will override this method
   *  @return a RefineResult enum such that:
   *          - REFINE_NONE if refinement was not needed, e.g. found a concrete
   * CEX
   *          - REFINE_SUCCESS if successfully refined / ruled out abstract CEX
   *          - REFINE_FAIL if failed during refinement
   *  NOTE the counterexample trace is accessible through cex_ which is
   *  set by block_all when a trace is found
   */
  virtual RefineResult refine()
  {
    // by default no refinement done -- e.g. assuming counterexamples are always
    // concrete
    return REFINE_NONE;
  }

  /** Check if a transition from the frontier can result in a bad state
   *  @param out an IC3Formula to populate with a state in the frontier
   *         that can reach bad in one step
   *  @return true iff bad is reachable from a state in the frontier
   */
  virtual bool reaches_bad(IC3Formula & out);

  // ********************************** Common Methods
  // These methods are common to all flavors of IC3 currently implemented

  /** Perform a IC3 step
   *  @param i
   */
  ProverResult step(int i);

  /** Perform the base IC3 step (zero case)
   */
  ProverResult step_01();

  /** Do a relative inductiveness check at frame i-1
   *  aka see if c at frame i is reachable from frame i-1
   *  @requires c -> F[i]
   *  @param i the frame number
   *  @param c the IC3Formula to check
   *  @param out the output collateral:
   *         if query is UNSAT and ic3_unsatcore_gen is true it will
   *           do a cheap unsat-core based generalization of c and set
   *           out to a subset (as a conjunction still)
   *           NOTE: depending on the IC3 variant and particular problem
   *           this can either be helpful or over-generalize
   *           it's worth experimenting with this option
   *         if query is SAT and get_pred is TRUE, will set out to
   *           the predecessor CTI after generalizing the predecessor
   *           (if that option is enabled)
   *         NOTE: this code does not call the (typically more expensive)
   *               inductive_generalization
   *  @param get_pred if set to true, will set out to the predecessor on a
   *                  SAT query. This option exists because rel_ind_check
   *                  is used in several places including the default
   *                  inductive_generalization procedure, where we don't need
   *                  predecessors if it's SAT. (in that case, SAT just means
   *                  we can't drop the literal we tried dropping)
   *  @return true iff c is inductive relative to frame i-1
   *  @ensures returns false  : out -> F[i-1] /\ \forall s in out . (s, c) \in
   * [T] returns true   : out unchanged, F[i-1] /\ T /\ c' is unsat
   */
  bool rel_ind_check(size_t i,
                     const IC3Formula & c,
                     IC3Formula & out,
                     bool get_pred = true);

  // Helper methods

  /** Attempt to block all proof goals
   *  to ensure termination, always choose proof goal with
   *  smallest time
   *  @return true iff all proof goals were blocked
   *  if returns false, sets cex_ to the trace
   */
  bool block_all();

  /** Check if the given proof goal is already blocked
   *  @param pg the proof goal
   *  @return true iff the proof goal is already blocked
   */
  bool is_blocked(const ProofGoal * pg);

  /** Try propagating all clauses from frame index i to the next frame.
   *  @param i the frame index to propagate
   *  @return true iff all the clauses are propagated (this means property was
   * proven)
   */
  bool propagate(size_t i);

  /** Calls predecessor_generalization to generalize the current
   *  model (assumes the current context is satisfiable)
   *  Then if approx_pregen_ is true will do a solver call
   *  to fix the predecessor if it was overgeneralized
   *  @requires rel_ind_check(i, c)
   *  @requires the solver_ context is currently satisfiable
   *  @param i the frame number
   *  @param c the term that failed to be blocked at i
   *         (over current state vars)
   *  @param pred the predecessor. originally passed as full assignment
   *         and then is updated in to be a generalized predecessor
   *  @ensures no calls to the solver_ because the context is polluted with
   *           other assertions
   */
  void predecessor_generalization_and_fix(size_t i,
                                          const smt::Term & c,
                                          IC3Formula & pred);

  /** Add a new frame
   *  If not F[0] then starts the frame with the property (implicitly)
   *  by assuming frame_label -> prop
   */
  void push_frame();

  /** Adds a constraint to frame i and (implicitly) all frames below it
   *  @param i highest frame to add constraint to
   *  @param constraint the constraint to add
   *  @param new_contraint true iff the constraint is a
   *         newly learned blocking constraint. In true, then subsumption check
   *         is performed
   */
  void constrain_frame(size_t i,
                       const IC3Formula & constraint,
                       bool new_constraint = true);

  /** Adds an implication frame_label_[i] -> constraint
   *  used as a helper in constrain_frame and when resetting solver
   *  to re-add those assertions
   *  @param i highest frame to add constraint to
   *  @param constraint the constraint associate with frame_label_[i]
   */
  void constrain_frame_label(size_t i, const IC3Formula & constraint);

  /** Add all the terms at Frame i
   *  Note: the frames_ data structure keeps terms only in the
   *  highest frame where they are known to hold
   *  Thus, asserting Fi actually needs to add terms
   *  from Fi and all frames after it
   *  @param i the frame number
   */
  void assert_frame_labels(size_t i) const;

  smt::Term get_frame_term(size_t i) const;

  void assert_trans_label() const;

  /** Check if there are common assignments
   *  between A and B
   *  i.e. if A /\ B is SAT
   *  @param A the first term
   *  @param B the second term
   *  @return true iff there is an intersection
   */
  bool check_intersects(const smt::Term & A, const smt::Term & B);

  /** Check if the term intersects with the initial states
   *  syntactic sugar for intersects(ts_.init(), t);
   *  @param t the term to check
   *  @return true iff t intersects with the initial states
   */
  bool check_intersects_initial(const smt::Term & t);

  void fix_if_intersects_initial(smt::TermVec & to_keep,
                                 const smt::TermVec & rem);

  /** Returns the highest frame this unit can be pushed to
   *  @param i the starting frame index
   *  @param u the IC3Formula to check how far it can be pushed
   *         automatically gets generalized with unsat-cores
   *  @return index >= i such that this unit can be added
   *          to that frame
   */
  size_t find_highest_frame(size_t i, IC3Formula & u);

  /** Returns a vector of equalities between input variables
   *  and their values in the current model
   *  @require solver_ state to be SAT
   *  @return vector of equalities
   */
  smt::TermVec get_input_values() const;

  /** Returns a vector of equalities between next state variables
   *  and their values in the current model
   *  @require solver_ state to be SAT
   *  @return vector of equalities
   */
  smt::TermVec get_next_state_values() const;

  /** Constructs a counterexample trace by following the
   *  next pointers of pg until the end of the trace
   *  @param pg the proof goal that starts at the initial states
   *  @param out the vector to populate
   *  NOTE: because of the reaches_bad implementation, the last goal
   *  in the pg chain is not a bad state, but rather a state that
   *  can reach bad in one step.
   *  thus this method always adds bad_ to the end of the vector
   */
  void reconstruct_trace(const ProofGoal * pg, smt::TermVec & out);

  /** Creates a reduce and of the vector of boolean terms
   *  It also sorts the vector by the hash
   *  Note: this will fail for empty vectors
   *  @param vec the vector of boolean terms
   *  @param slv (optional) the solver to use, defaults to solver_
   *  @return the conjunction of all the terms
   */
  smt::Term make_and(smt::TermVec vec, smt::SmtSolver slv = NULL) const;

  /** Pushes a solver context and keeps track of the context level
   *  updates solver_context_
   */
  inline void push_solver_context()
  {
    solver_->push();
    solver_context_++;
  }

  /** Pops a solver context and keeps track of the context level
   *  updates solver_context_
   */
  inline void pop_solver_context()
  {
    solver_->pop();
    solver_context_--;
  }

  inline smt::Result check_sat()
  {
    num_check_sat_since_reset_++;
    return solver_->check_sat();
  }

  inline smt::Result check_sat_assuming(const smt::TermVec & assumps)
  {
    num_check_sat_since_reset_++;
    return solver_->check_sat_assuming(assumps);
  }

  /** Attempts to reset the solver and re-add constraints
   *  NOTE: not all solvers support reset_assertions, in which case the
   * exception is just caught and things continue on as normal
   */
  virtual void reset_solver();

  inline size_t frontier_idx() const { return frames_.size() - 1; }

  /** Create a boolean label for a given term
   *  These are cached in labels_
   *  good for using unsat cores
   *
   *  @param t a boolean formula to create a label for
   *  @return the indicator variable label for this term
   */
  smt::Term label(const smt::Term & t);

  /** Returns true iff this label
   *  already has semantics assigned at the base context level
   */
  virtual bool is_global_label(const smt::Term & l) const;

  /** Negates a term by stripping the leading Not if it's there,
   ** or applying Not if the term is not already negated.
   */
  smt::Term smart_not(const smt::Term & t) const;
};

}  // namespace pono
