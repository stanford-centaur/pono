/*********************                                                        */
/*! \file interpolantmc.h
 ** \verbatim
 ** Top contributors (to current version):
 **   Makai Mann, Ahmed Irfan
 ** This file is part of the pono project.
 ** Copyright (c) 2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief A straightforward implementation of interpolation based model
 *checking.
 **        See Interpolation and SAT-based Model Checking
 **
 **/

#pragma once

#include <cstdint>

#include "engines/prover.h"
#include "options/options.h"
#include "smt-switch/smt.h"

namespace pono {

class InterpolantMC : public SafetyProver
{
 public:
  InterpolantMC(const SafetyProperty & p,
                const TransitionSystem & ts,
                const smt::SmtSolver & solver,
                PonoOptions opt = {},
                Engine engine = Engine::INTERP);

  typedef SafetyProver super;

  void initialize() override;

  void reset_env() override;

  ProverResult check_until(int k) override;

 protected:
  bool step(const int i);
  bool step_0();

  /**
   * @brief Check whether the reached states have converged to a fixed point,
   * that is, whether the newly computed interpolant is already covered by the
   * reached states.
   *
   * This method modifies the solver stack.
   *
   * @param new_itp the newly computed interpolant
   * @param reached The reached states, represented as the disjunction of
   *        previously computed interpolants.
   * @param interp_count The number of interpolants computed at the current
   *        unrolling step.
   * @return true iff `new_itp` is already covered by `reached`
   */
  bool has_converged(const smt::Term & new_itp,
                     const smt::Term & reached,
                     const int & interp_count);

  /**
   * @brief Make the step of the suffix from time k: trans@k if the transition
   * relation is right-total, and otherwise bad@k \/ trans@k, as a path to a
   * bad state needs a step only before it reaches the bad state.
   *
   * Unless the transition relation is right-total, a bad state need not have
   * a successor. With the latter steps, transB_ /\ bad@i holds from exactly
   * the states that reach a bad state within the suffix up to time i, like
   * bad@1 \/ (trans@1 /\ (bad@2 \/ (trans@2 /\ ...))), as a path that
   * reaches a bad state earlier can stay in it. So every interpolant excludes
   * the bad states, and the fixed point implies the property.
   *
   * @param k the time of the state that the step leaves
   * @return the step from time k to time k+1
   */
  smt::Term suffix_step(int k);

  // configurable options
  const bool use_frontier_simpl_;
  const InterpPropsEnum interp_props_;
  const bool unroll_eagerly_;
  const bool interp_backward_;

  smt::SmtSolver interpolator_;
  // for translating terms to interpolator_
  smt::TermTranslator to_interpolator_;
  // for translating terms to solver_
  smt::TermTranslator to_solver_;

  // set to true when a concrete_cex is found
  bool concrete_cex_ = false;

  bool trans_total_ = false;  ///< whether ts_ is known to be right-total

  smt::Term init0_;
  smt::Term transA_;
  smt::Term transB_;  ///< the steps of the suffix, see suffix_step()
  /// a disjunction of bads in the suffix, or only its last bad if the
  /// transition relation is not right-total, see suffix_step()
  smt::Term bad_disjuncts_;

#ifndef NDEBUG
  std::uint32_t total_interp_call_count_ = 0;
  double total_interp_call_time_ = 0.0;
#endif

};  // class InterpolantMC

}  // namespace pono
