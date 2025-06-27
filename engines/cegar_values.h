/*********************                                                        */
/*! \file cegar_values.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief A simple CEGAR loop that abstracts values with frozen variables
**        and refines by constraining the variable to the value again
**
**/

#pragma once

#include "core/unroller.h"
#include "engines/cegar.h"

namespace pono {
template <class Prover_T>
class CegarValues : public CEGAR<Prover_T>
{
  typedef CEGAR<Prover_T> super;

 public:
  CegarValues(const SafetyProperty & p,
              const TransitionSystem & ts,
              const smt::SmtSolver & solver,
              PonoOptions opt = PonoOptions());

  ProverResult check_until(int k) override;

  void initialize() override;

 protected:
  TransitionSystem conc_ts_;
  TransitionSystem & prover_ts_;

  // solver and associated infrastructure for
  // unrolling based refinement
  smt::SmtSolver cegval_solver_;
  smt::TermTranslator to_cegval_solver_;
  smt::TermTranslator from_cegval_solver_;
  TransitionSystem cegval_ts_;
  Unroller cegval_un_;
  smt::UnorderedTermMap to_vals_;
  smt::Term cegval_bad_;

  smt::UnorderedTermMap cegval_labels_;  // labels for each abstract value

  void cegar_abstract() override;

  bool cegar_refine() override;

  void refine_subprover_ts(const smt::UnorderedTermSet & axioms);
};

}  // namespace pono
