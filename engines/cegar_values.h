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

#include "engines/cegar.h"

namespace pono {
template <class Prover_T>
class CegarValues : public CEGAR<Prover_T>
{
  typedef CEGAR<Prover_T> super;

 public:
  CegarValues(const Property & p,
              const TransitionSystem & ts,
              const smt::SmtSolver & solver,
              PonoOptions opt = PonoOptions());

  ProverResult check_until(int k) override;

  void initialize() override;

 protected:
  TransitionSystem conc_ts_;
  TransitionSystem & prover_ts_;
  smt::UnorderedTermMap to_vals_;

  void cegar_abstract() override;

  bool cegar_refine() override;
};

}  // namespace pono
