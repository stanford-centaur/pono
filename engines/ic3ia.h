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

#include "core/unroller.h"
#include "engines/ic3.h"
#include "modifiers/implicit_predicate_abstractor.h"
#include "smt-switch/term_translator.h"

namespace pono {

class IC3IA : public IC3
{
 public:
  // itp_se is the SolverEnum for the interpolator

  IC3IA(const Property & p,
        const TransitionSystem & ts,
        const smt::SmtSolver & s,
        PonoOptions opt = PonoOptions());

  virtual ~IC3IA() {}

  typedef IC3 super;

 protected:
  // Note: important that conc_ts_ and abs_ts_ are before ia_
  //       because we will pass them to ia_ and they must be
  //       be initialized first

  TransitionSystem conc_ts_;

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

  // Since MathSAT is the best solver for IC3IA it helps to use
  // its bool_model_generation option which doesn't enable
  // model generation for theories
  // instead, we assign indicator labels to check the value of
  // predicates. This should still work fine with other solvers
  smt::UnorderedTermMap lbl2pred_;
  smt::UnorderedTermSet predlbls_;
  smt::UnorderedTermSet all_lbls_;  ///< for debugging assertions only
                                    ///< keeps track of both polarities
                                    ///< mostly needed because some solvers
                                    ///< do automatic top-level propagation
                                    ///< which means you can't count on a symbol
                                    ///< staying a symbol

  // Hacked in for CVC4 SyGuS predicate experimentation

  // need to be able to unroll abstract ts (regular ic3ia doesn't)
  // and currently the unroller_ is over the conc_ts_
  Unroller abs_unroller_;

  smt::UnorderedTermSet max_terms_;  ///< largest non-Boolean terms in TS

  // extra members for this hacked in stuff
  std::unordered_set<smt::SortKind>
      all_sort_kinds_;  ///< all sort kinds appearing in TS

  smt::TermVec pred_candidates_;  ///< known predicates not in predset_
                                  ///< these might be able to rule out
                                  ///< abstract counterexamples

  /** Overriding the method. This will return the concrete_ts_ because ts_ is an
   *  abstraction of concrete_ts_.
   */
  TransitionSystem & prover_interface_ts() override { return conc_ts_; };
  // pure virtual method implementations

  IC3Formula get_model_ic3formula() const override;

  bool ic3formula_check_valid(const IC3Formula & u) const override;

  // need to override this because IC3IA is not as restricted as
  // (bit-level) IC3
  void check_ts() const override;

  void initialize() override;

  void abstract() override;

  RefineResult refine() override;

  void reset_solver() override;

  bool is_global_label(const smt::Term & l) const override;

  // specific to IC3IA

  void reabstract();

  /** Adds predicate to abstraction
   *  (calls ia_.add_predicate)
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

  /** Synthesize predicates using CVC4 SyGuS
   *  used as a helper function for cvc4_find_preds
   *  @param abs_trace the unrolled abstract trace (over solver_ terms)
   *  @param state variables over solver_ terms
   *         will respect this order of state variables (that's why we can't
   *         just get the set of state variables from the TS)
   *  @param unrolled_var_args - vector of pairs where first is unrolled
   *         next vars and second is unrolled abstract variables
   *         (over solver_ terms)
   *  @param free_vars - set of all free variables in abs_trace (over solver_
   *         terms). Includes unrolled input variables also.
   *  @param num_preds - how many predicates to look for
   *  @param out_preds - set to add synthesized predicates to
   *  @return true iff predicates were found that rule out this abstract trace
   */
  bool cvc4_synthesize_preds(
      const smt::Term & abs_trace,
      const smt::TermVec & statevars,
      const std::vector<std::pair<smt::TermVec, smt::TermVec>> &
          unrolled_var_args,
      const smt::UnorderedTermSet & free_vars,
      size_t num_preds,
      smt::UnorderedTermSet & out_preds);
};

}  // namespace pono
