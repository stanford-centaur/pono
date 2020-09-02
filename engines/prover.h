/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Makai Mann, Ahmed Irfan, Florian Lonsing
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

#include "core/prop.h"
#include "core/proverresult.h"
#include "core/ts.h"
#include "core/unroller.h"
#include "options/options.h"

#include "smt-switch/smt.h"

namespace pono {
class Prover
{
 public:
  Prover(const Property & p, smt::SolverEnum se);
  Prover(const Property & p, const smt::SmtSolver & s);
  Prover(const PonoOptions & opt, const Property & p, smt::SolverEnum se);
  Prover(const PonoOptions & opt, const Property & p, const smt::SmtSolver & s);

  virtual ~Prover();

  virtual void initialize();

  virtual ProverResult check_until(int k) = 0;

  bool witness(std::vector<smt::UnorderedTermMap> & out);

  virtual ProverResult prove();

 protected:
  smt::SmtSolver solver_;
  smt::TermTranslator to_prover_solver_;
  const Property property_;
  const TransitionSystem &
      ts_;  ///< convenient reference to transition system in property
  const TransitionSystem &
      orig_ts_;  ///< reference to original TS before copied to new solver

  Unroller unroller_;

  int reached_k_;

  smt::Term bad_;

  PonoOptions options_;
};
}  // namespace pono
