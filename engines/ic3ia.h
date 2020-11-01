/*********************                                                  */
/*! \file ic3ia.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief IC3 via Implicit Predicate Abstraction (IC3IA) implementation
**        based on
**
**        IC3 Modulo Theories via Implicit Predicate Abstraction
**            -- Alessandro Cimatti, Alberto Griggio,
**               Sergio Mover, Stefano Tonetta
**
**        and the open source implementation:
**
**        https://es-static.fbk.eu/people/griggio/ic3ia/index.html
**/

#include "engines/mbic3.h"
#include "modifiers/implicit_predicate_abstractor.h"

namespace pono {

// process is mostly the same as model based IC3,
// but overrides a few methods to only use cubes/clauses over predicates
class IC3IA : public ModelBasedIC3
{
 public:
  IC3IA(Property & p, smt::SolverEnum se);
  IC3IA(Property & p, const smt::SmtSolver & slv);
  IC3IA(const PonoOptions & opt, Property & p, smt::SolverEnum se);
  IC3IA(const PonoOptions & opt, Property & p, const smt::SmtSolver & slv);
  virtual ~IC3IA();

  typedef ModelBasedIC3 super;

 protected:
  void initialize() override;

  bool block_all() override;

  // modified to get a bad cube over the predicates
  bool intersects_bad() override;

  // modified to get a bad cube over the predicates
  // TODO: refactor so that generalization and getting cube are different
  bool get_predecessor(size_t i,
                       const Conjunction & c,
                       Conjunction & out_pred) override;

  // use to get a predicate assignment from the current solver context
  Conjunction get_conjunction_from_model();

  /** Overrides the set labels method
   *  to instead use the abstract transition relation
   */
  void set_labels() override;

  /** Adds predicate to abstraction
   *  (calls ia_.add_predicate)
   *  and also incrementally updates the local transition relation
   *  and declares a new predicate state var (in pred_statevars_)
   *  @param pred the predicate over current state variables
   */
  void add_predicate(const smt::Term & pred);

  /** Refines an abstract counterexample trace
   *  by looking for new predicates
   *  @param pg the proof goal that reached the initial state
   *  @return FALSE if there's really a counterexample, TRUE if
   *          the counterexample was refined and UNKNOWN if no
   *          new predicates were found
   */
  ProverResult refine(ProofGoal pg);

  // TODO figure out relevant state

  RelationalTransitionSystem abs_ts_;

  ImplicitPredicateAbstractor ia_;

  smt::TermVec pred_statevars_;
  ///< state variables that IC3IA sees
  ///< these represent predicates and should stay aligned with the
  ///< predicate vector in the abstraction (ia_.predicates())
  ///< e.g. pred_statevars_[i] corresponds to ia_.predicates()[i]

  // TODO: remove if not using them
  // smt::UnorderedTermMap predstate2next_;
  // ///< maps predicate state variable to the next state version

  // smt::UnorderedTermMap nextpredstate2curr_;
  // ///< maps predicate next state variable to the current version

  // useful sorts
  smt::Sort boolsort_;
};

}  // namespace pono
