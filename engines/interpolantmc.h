/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Makai Mann, Ahmed Irfan
 ** This file is part of the pono project.
 ** Copyright (c) 2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief
 **
 **
 **/

#pragma once

#include "engines/prover.h"

#include "smt-switch/smt.h"

namespace pono {

class InterpolantMC : public Prover
{
 public:
  // IMPORTANT: assume the property was built using the interpolating solver
  InterpolantMC(const Property & p, smt::SmtSolver & slv, smt::SmtSolver & itp);
  InterpolantMC(const PonoOptions & opt,
                const Property & p,
                smt::SmtSolver & slv,
                smt::SmtSolver & itp);
  ~InterpolantMC();

  typedef Prover super;

  void initialize() override;

  ProverResult check_until(int k) override;

 private:
  bool step(int i);
  bool step_0();

  void reset_assertions(smt::SmtSolver & s);

  bool check_entail(const smt::Term & p, const smt::Term & q);

  smt::SmtSolver & interpolator_;
  // for translating terms to interpolator_
  smt::TermTranslator to_interpolator_;
  // for translating terms to solver_
  smt::TermTranslator to_solver_;

  // set to true when a concrete_cex is found
  bool concrete_cex_;

  smt::Term init0_;
  smt::Term transA_;
  smt::Term transB_;

};  // class InterpolantMC

}  // namespace pono
