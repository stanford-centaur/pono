/*********************                                                  */
/*! \file interp_seq_mc.h
 ** \verbatim
 ** Top contributors (to current version):
 **   Po-Chun Chien
 ** This file is part of the pono project.
 ** Copyright (c) 2019 by the authors listed in the file AUTHORS
 ** (in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief Implementation of Dual Approximated Reachability.
 **
 ** This unbounded model-checking algorithm was introduced by Yakir Vizel,
 ** Orna Grumberg, and Sharon Shoham in the paper "Intertwined Forward-Backward
 ** Reachability Analysis Using Interpolants" published in TACAS 2013
 ** (https://doi.org/10.1007/978-3-642-36742-7_22).
 **
 **/

#pragma once

#include "engines/prover.h"
#include "smt-switch/smt.h"

namespace pono {

class DualApproxReach : public SafetyProver
{
 public:
  DualApproxReach(const SafetyProperty & p,
                  const TransitionSystem & ts,
                  const smt::SmtSolver & slv,
                  PonoOptions opt = PonoOptions());

  ~DualApproxReach();

  typedef SafetyProver super;

  void initialize() override;

  void reset_env() override;

  ProverResult check_until(int k) override;

 protected:
  bool step(int i);
  bool step_0();

  bool local_strengthen();
  bool global_strengthen();
  void pairwise_strengthen(const size_t idx);

  void update_term_map(size_t i);

  bool check_fixed_point();
  bool check_fixed_point(const smt::TermVec & reach_seq,
                         smt::Term & fixed_point);
  bool check_entail(const smt::Term & p, const smt::Term & q);

  smt::SmtSolver interpolator_;
  // for translating terms to interpolator_
  smt::TermTranslator to_interpolator_;
  // for translating terms to solver_
  smt::TermTranslator to_solver_;

  // set to true when a concrete_cex is found
  bool concrete_cex_;

  // forward reachability sequence: <F_0=Init, F_0, F_1, ...>
  smt::TermVec forward_seq_;
  // backward reachability sequence: <B_0=Bad, B_1, B_2, ...>
  smt::TermVec backward_seq_;
};  // class DualApproxReach

}  // namespace pono
