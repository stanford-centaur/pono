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
#include "smt-switch/smt.h"

namespace pono {

class InterpolantMC : public SafetyProver
{
 public:
  InterpolantMC(const SafetyProperty & p,
                const TransitionSystem & ts,
                const smt::SmtSolver & slv,
                PonoOptions opt = PonoOptions());

  ~InterpolantMC();

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

  smt::Term init0_;
  smt::Term transA_;
  smt::Term transB_;
  smt::Term bad_disjuncts_;  ///< a disjunction of bads in the suffix

#ifndef NDEBUG
  std::uint32_t total_interp_call_count_ = 0;
  double total_interp_call_time_ = 0.0;
#endif

};  // class InterpolantMC

}  // namespace pono
