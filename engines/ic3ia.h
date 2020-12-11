/*********************                                                  */
/*! \file ic3ia.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan
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
**
**  within Pono, we are building on the bit-level IC3 instead of directly
**  on IC3Base, because a lot of the functionality is the same
**  In particular, we don't need to override either of the generalization
**  functions. Instead focusing on abstract/refine.
**
**/

#pragma once

#include "engines/ic3.h"
#include "modifiers/implicit_predicate_abstractor.h"
#include "smt-switch/term_translator.h"

namespace pono {

class IC3IA : public IC3
{
 public:
  // itp_se is the SolverEnum for the interpolator

  IC3IA(Property & p,
        smt::SolverEnum se,
        smt::SolverEnum itp_se = smt::SolverEnum::MSAT);
  IC3IA(Property & p,
        const smt::SmtSolver & s,
        smt::SolverEnum itp_se = smt::SolverEnum::MSAT);
  IC3IA(const PonoOptions & opt,
        Property & p,
        smt::SolverEnum se,
        smt::SolverEnum itp_se = smt::SolverEnum::MSAT);
  IC3IA(const PonoOptions & opt,
        Property & p,
        const smt::SmtSolver & s,
        smt::SolverEnum itp_se = smt::SolverEnum::MSAT);
  virtual ~IC3IA() {}

  typedef IC3 super;

 protected:
  // Note: important that conc_ts_ and abs_ts_ are before ia_
  //       because we will pass them to ia_ and they must be
  //       be initialized first

  TransitionSystem & conc_ts_;  ///< convenient reference to the concrete ts

  RelationalTransitionSystem abs_ts_;  ///< the abstract ts
                                       ///< after initialize, ts_ will point to
                                       ///< this because the methods from IC3
                                       ///< should operate on the abstraction

  ImplicitPredicateAbstractor ia_;

  smt::UnorderedTermSet predset_;  ///< set of current predicates
  // useful for checking if predicate has been added already
  // also available as a vector in ia_.predicates()

  smt::SmtSolver interpolator_;  ///< interpolator for refinement
  smt::TermTranslator
      to_interpolator_;  ///< transfer terms from solver_ to interpolator_
  smt::TermTranslator
      to_solver_;  ///< transfer terms from interpolator_ to solver_

  size_t longest_cex_length_;  ///< keeps track of longest (abstract)
                               ///< counterexample

  // HACK
  // hacked in to evaluate CVC4
  // if done for real, should be sure to do this OR the interpolator, not both
  smt::SmtSolver cvc4_;
  smt::TermTranslator to_cvc4_;
  smt::TermTranslator from_cvc4_;

  // pure virtual method implementations

  IC3Formula get_ic3_formula(smt::TermVec * inputs = nullptr,
                             smt::TermVec * nexts = nullptr) const override;

  bool ic3_formula_check_valid(const IC3Formula & u) const override;

  // need to override this because IC3IA is not as restricted as
  // (bit-level) IC3
  void check_ts() const override;

  void initialize() override;

  void abstract() override;

  RefineResult refine() override;

  // specific to IC3IA

  /** Adds predicate to abstraction
   *  (calls ia_.add_predicate)
   *  and also incrementally updates the local transition relation
   *  and declares a new predicate state var (in pred_statevars_)
   *  @param pred the predicate over current state variables
   *  @return true iff the predicate was new (not seen before)
   */
  bool add_predicate(const smt::Term & pred);

  /** Register a state variable mapping in to_solver_
   *  This is a bit ugly but it's needed because symbols aren't created in
   * to_solver_ so it needs the mapping from interpolator_ symbols to solver_
   * symbols
   *  TODO look into a cleaner solution
   *  @param i the unrolling for state variables
   *         makes sure not to repeat work
   */
  void register_symbol_mappings(size_t i);

  // Hacked in to experiment with CVC4

  /** Given a counterexample trace (over state vars)
   *  Unroll the trace and ask CVC4 SyGuS for predicate(s)
   *  That makes the abstract trace unsat
   *  @param cex a vector storing the state variable assignments at each step of
   * an abstract trace
   *  @param out_preds the set to add predicates to
   *  @param return true if a predicate was found
   */
  bool cvc4_find_preds(const smt::TermVec & cex,
                       smt::UnorderedTermSet & out_preds);
};

}  // namespace pono
