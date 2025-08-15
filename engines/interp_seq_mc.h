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
 ** \brief Implementation of interpolation-sequence based model checking.
 **
 ** This unbounded model checking algorithm was introduced by Yakir Vizel and
 ** Orna Grumberg in FMCAD 2009 (https://doi.org/10.1109/FMCAD.2009.5351148).
 **
 **/

#pragma once

#include "engines/prover.h"
#include "smt-switch/smt.h"

namespace pono {

class InterpSeqMC : public SafetyProver
{
 public:
  InterpSeqMC(const SafetyProperty & p,
              const TransitionSystem & ts,
              const smt::SmtSolver & slv,
              PonoOptions opt = PonoOptions());

  ~InterpSeqMC();

  typedef SafetyProver super;

  void initialize() override;

  void reset_env() override;

  ProverResult check_until(int k) override;

 protected:
  bool step(int i);
  bool step_0();

  void update_term_map(size_t i);

  bool check_fixed_point();
  bool check_entail(const smt::Term & p, const smt::Term & q);
  void check_itp_sequence(const smt::TermVec & int_formulas,
                          const smt::TermVec & int_itp_seq);

  smt::SmtSolver interpolator_;
  // for translating terms to interpolator_
  smt::TermTranslator to_interpolator_;
  // for translating terms to solver_
  smt::TermTranslator to_solver_;

  // set to true when a concrete_cex is found
  bool concrete_cex_;

  // reachability sequence: <R_0=Init, R_1, R_2, ...>
  smt::TermVec reach_seq_;
  // transition at each time step: <Init(0) & TR(0, 1), TR(1, 2), ...>
  // note that at 0th step, Init is conjoined with TR
  smt::TermVec trans_seq_;
  smt::TermVec int_trans_seq_;

};  // class InterpSeqMC

}  // namespace pono
