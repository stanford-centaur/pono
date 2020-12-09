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
**           - implement get_ic3_formula and give it semantics to produce the
**             corresponding IC3Formula for your flavor of IC3
**             (assumes solver_'s state is SAT from a failed rel_ind_check)
**             also need to be able to give model (as formulas) for input values
**             and next-state variable values if inputs/nexts are non-null,
*respectively
**           - implement inductive_generalization
**           - implement generalize_predecessor
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
**           - implement all the ic3_formula_* functions for creating and
**             manipulating an IC3Formula if the defaults are not right
**           - implement abstract() and refine() if this is a CEGAR
**             flavor of IC3
**
**        Important Notes:
**           - be sure to use [push/pop]_solver_context instead of using
**             the solver interface directly so that the assertions on
**             the solver context (tracked externally) are correct
**
**/
#pragma once

#include "engines/prover.h"
#include "smt-switch/utils.h"

namespace pono {

struct IC3Formula
{
  // nullary constructor
  IC3Formula() : term(nullptr) {}
  IC3Formula(const smt::Term & t, const smt::TermVec & c, bool n)
      : term(t), children(c), disjunction(n)
  {
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

  bool is_disjunction() const { return disjunction; };

  smt::Term term;  ///< term representation of this formula
  smt::TermVec
      children;  ///< flattened children of either a disjunction or conjunction
  bool disjunction;  ///< true if currently representing a disjunction
  // NOTE: treating the disjunction as the primary orientation
  //       e.g. aligned with what's kept in frames (disjunctions), not proof
  //       goals (conjunctions), but IC3Formula can represent both
};

// TODO change back to ProofGoal once refactor is done
// don't want to clash with name in MBIC3 for now
struct IC3Goal
{
  // based on open-source ic3ia ProofObligation
  IC3Formula target;
  size_t idx;
  // TODO: see if we can make this a unique_ptr
  //       made it complicated to move from this struct to another place
  std::shared_ptr<IC3Goal> next;

  // null constructor
  IC3Goal() : idx(0), next(nullptr) {}

  IC3Goal(IC3Formula u, size_t i, const std::shared_ptr<IC3Goal> & n)
      : target(u), idx(i), next(n)
  {
  }

  IC3Goal(const IC3Goal & other)
      : target(other.target), idx(other.idx), next(other.next)
  {
  }
};

class IC3Base : public Prover
{
 public:
  /** IC3Base constructors take the normal arguments for a prover
   *  + a function that can create an IC3Formula
   *  Depending on the derived class IC3 implementation, the exact
   *  type of IC3Formula will differ: e.g. Clause, Disjunction
   */
  IC3Base(Property & p, smt::SolverEnum se);
  IC3Base(Property & p, const smt::SmtSolver & s);
  IC3Base(const PonoOptions & opt, Property & p, smt::SolverEnum se);
  IC3Base(const PonoOptions & opt, Property & p, const smt::SmtSolver & s);

  typedef Prover super;

  void initialize() override;

  ProverResult check_until(int k) override;

  bool witness(std::vector<smt::UnorderedTermMap> & out) override;

 protected:

  smt::UnsatCoreReducer reducer_;

  ///< keeps track of the current context-level of the solver
  // NOTE: if solver is passed in, it could be off
  //       currently no way to check
  //       then solver_context_ is relative to the starting context
  size_t solver_context_;

  ///< the frames data structure.
  ///< a vector of the given Unit template
  ///< which changes depending on the implementation
  std::vector<std::vector<IC3Formula>> frames_;

  ///< stack of outstanding proof goals
  std::vector<IC3Goal> proof_goals_;

  // labels for activating assertions
  smt::Term init_label_;       ///< label to activate init
  smt::Term trans_label_;      ///< label to activate trans
  smt::TermVec frame_labels_;  ///< labels to activate frames
  smt::UnorderedTermMap labels_;  //< labels for unsat cores

  IC3Goal cex_pg_;  ///< if a proof goal is traced back to init
                    ///< this gets set to the first proof goal
                    ///< in the trace
                    ///< otherwise starts null, can check that
                    ///< cex_pg_.target.term is a nullptr

  // useful terms
  smt::Term solver_true_;

  // TODO Make sure all comments are updated!

  // *************************** Main Methods *********************************

  // ************************** Virtual Methods *******************************
  // IMPORTANT for derived classes
  // These methods should be implemented by a derived class for a particular
  // "flavor" of IC

  // Pure virtual methods that must be overridden

  /** Get an IC3Formula from the current model
   *  @requires last call to check_sat of solver_ was satisfiable and context
   * hasn't changed
   *  @param inputs - pointer to a vector. If non-null populate with input
   *                  variable model
   *  @param nexts - pointer to a vector. If non-null populate with next
   *                 state variable model
   *  @return an IC3Formula over current state variables with is_disjunction
   *          false depending on the flavor of IC3, this might be a boolean
   * cube, a theory cube, a cube of predicates, etc... AND if inputs non-null,
   * then include model for inputs, e.g. as equalities if nexts non-null, then
   * include model for next state vars, e.g. add equalities to vector
   */
  virtual IC3Formula get_ic3_formula(smt::TermVec * inputs = nullptr,
                                     smt::TermVec * nexts = nullptr) const = 0;

  /** Check whether a given IC3Formula is valid
   *  e.g. if this is a boolean clause it would
   *    check that it's a disjunction of literals
   *  (for debugging)
   *  @param u the IC3Formula to check
   *  @return true iff this is a valid IC3Formula for this
   *          flavor of IC3
   */
  virtual bool ic3_formula_check_valid(const IC3Formula & u) const = 0;

  /** Attempt to generalize before blocking in frame i
   *  The standard approach is inductive generalization
   *  @requires !rel_ind_check(i, c, _)
   *  @param i the frame number to generalize it against
   *  @param c the IC3Formula that should be blocked
   *  @return a vector of IC3Formulas interpreted as a conjunction of
   * IC3Formulas. Standard IC3 implementations will have a size one vector (e.g.
   * a single clause) Let the returned conjunction term be d
   *  @ensures d -> !c and F[i-1] /\ d /\ T /\ !d' is unsat
   *           e.g. it blocks c and is inductive relative to F[i-1]
   */
  virtual std::vector<IC3Formula> inductive_generalization(
      size_t i, const IC3Formula & c) = 0;

  /** Generalize a counterexample
   *  @requires rel_ind_check(i, c)
   *  @requires the solver_ context is currently satisfiable
   *  @param i the frame number
   *  @param c the IC3Formula to find a general predecessor for
   *  @return a new IC3Formula d
   *  @ensures d -> F[i-1] /\ forall s \in [d] exists s' \in [c]. (d,c) \in [T]
   *  @ensures no calls to the solver_ because the context is polluted with
   *           other assertions
   */
  virtual IC3Formula generalize_predecessor(size_t i, const IC3Formula & c) = 0;

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
  virtual IC3Formula ic3_formula_disjunction(const smt::TermVec & c) const;

  /** Creates a conjunction IC3Formula from a vector of terms
   *  @param c the children terms
   *  @ensures resulting IC3Formula children == c
   *  @ensures resulting IC3Formula with is_disjunction false
   *  e.g. for a ClauseHandler, this method will create a cube
   *  note: assumes the children are already in the right polarity
   *  (doesn't negate them)
   */
  virtual IC3Formula ic3_formula_conjunction(const smt::TermVec & c) const;

  /** Negates an IC3Formula
   *  @param u the IC3Formula to negate
   */
  virtual IC3Formula ic3_formula_negate(const IC3Formula & u) const;

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
    ;
  }

  /** Refines an abstract transition system
   *  By default, this does nothing (e.g. assumes ts_ is not abstract)
   *  Any CEGAR implementations of IC3 will override this method
   *  @return a RefineResult enum such that:
   *          - REFINE_NONE if refinement was not needed, e.g. found a concrete
   * CEX
   *          - REFINE_SUCCESS if successfully refined / ruled out abstract CEX
   *          - REFINE_FAIL if failed during refinement
   *  NOTE the counterexample trace is accessible through cex_pg_ which is
   *  set by block_all when a trace is found
   *  can reconstruct the trace (without input variable values) by following
   *  the IC3Goal next field iteratively
   */
  virtual RefineResult refine()
  {
    // by default no refinement done -- e.g. assuming counterexamples are always
    // concrete
    return REFINE_NONE;
  }

  /** Check if a transition from the second to last frame can result in a bad
   * state
   *  @return true iff the last frame intersects with bad
   *  post-condition: if true is returned, a bad IC3Formula is added to proof
   * goals This method can be overriden if you want to add more than a single
   *  IC3Formula that intersects bad to the proof goals
   */
  virtual bool intersects_bad();

  // ********************************** Common Methods
  // These methods are common to all flavors of IC3 currently implemented

  /** Perform a IC3 step
   *  @param i
   */
  ProverResult step(int i);

  /** Perform the base IC3 step (zero case)
   */
  ProverResult step_0();

  /** Do a relative inductiveness check at frame i-1
   *  aka see if c at frame i is reachable from frame i-1
   *  @requires c -> F[i]
   *  @param i the frame number
   *  @param c the IC3Formula to check
   *  @param out the output collateral: a vector interpreted as a conjunction of
   * IC3Formulas if the check succeeds (e.g. is UNSAT), then returns a vector of
   * blocking units to be added to Frame i if the check fails (e.g. is SAT),
   * then returns a vector of predecessors Note 1: this method calls
   * inductive_generalization and generalize_predecessor if options_.ic3_pregen_
   * and options_.ic3_indgen_ are set, respectively Note 2: in most cases, the
   * vector returned will be size one
   *  @return true iff c is inductive relative to frame i-1
   *  @ensures returns false  : out -> F[i-1] /\ \forall s in out . (s, c) \in
   * [T] returns true   : out unchanged, F[i-1] /\ T /\ c' is unsat
   */
  bool rel_ind_check(size_t i,
                     const IC3Formula & c,
                     std::vector<IC3Formula> & out);

  // Helper methods

  /** Attempt to block all proof goals
   *  to ensure termination, always choose proof goal with
   *  smallest time
   *  @return true iff all proof goals were blocked
   *  if returns false, sets cex_pg_ to the first IC3Goal
   *  of a trace, e.g. the trace can be recovered by following
   *  pg.next iteratively
   */
  bool block_all();

  /** Attempt to block cube c at frame i
   *  @param i the frame number
   *  @param c the cube to try blocking
   *  @return true iff the cube was blocked, otherwise a new proof goal was
   * added to the proof goals
   */
  bool block(const IC3Goal & pg);

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
  void constrain_frame(size_t i, const IC3Formula & constraint);

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

  /** Check if there are more proof goals
   *  @return true iff there are more proof goals
   */
  bool has_proof_goals() const;

  /** Gets a new proof goal (and removes it from proof_goals_)
   *  @requires has_proof_goals()
   *  @return a proof goal with the lowest available frame number
   *  @alters proof_goals_
   *  @ensures returned proof goal is from lowest frame in proof goals
   */
  IC3Goal get_next_proof_goal();

  /** Create and add a proof goal for cube c for frame i
   *  @param c the cube of the proof goal
   *  @param i the frame number for the proof goal
   *  @param n pointer to the proof goal that led to this one -- null for bad
   *  (i.e. end of trace)
   */
  void add_proof_goal(const IC3Formula & c,
                      size_t i,
                      std::shared_ptr<IC3Goal> n);

  /** Check if there are common assignments
   *  between A and B
   *  i.e. if A /\ B is SAT
   *  @param A the first term
   *  @param B the second term
   *  @return true iff there is an intersection
   */
  bool intersects(const smt::Term & A, const smt::Term & B);

  /** Check if the term intersects with the initial states
   *  syntactic sugar for intersects(ts_.init(), t);
   *  @param t the term to check
   *  @return true iff t intersects with the initial states
   */
  bool intersects_initial(const smt::Term & t);

  void fix_if_intersects_initial(smt::TermVec & to_keep,
                                 const smt::TermVec & rem);

  /** Returns the highest frame this unit can be pushed to
   *  @param i the starting frame index
   *  @param u the unit to check how far it can be pushed
   *  @return index >= i such that this unit can be added
   *          to that frame
   */
  size_t find_highest_frame(size_t i, const IC3Formula & u);

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
  void push_solver_context();

  /** Pops a solver context and keeps track of the context level
   *  updates solver_context_
   */
  void pop_solver_context();

  /** Create a boolean label for a given term
   *  These are cached in labels_
   *  good for using unsat cores
   *
   *  @param t a boolean formula to create a label for
   *  @return the indicator variable label for this term
   */
  smt::Term label(const smt::Term & t);

  /** Negates a term by stripping the leading Not if it's there,
   ** or applying Not if the term is not already negated.
   */
  smt::Term smart_not(const smt::Term & t) const;
};

}  // namespace pono
