/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Ahmed Irfan, Makai Mann
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

namespace pono {

class Bmc : public Prover
{
 public:
  Bmc(Property & p, smt::SolverEnum se);
  Bmc(Property & p, const smt::SmtSolver & solver);
  Bmc(const PonoOptions & opt, Property & p, smt::SolverEnum se);
  Bmc(const PonoOptions & opt,
      Property & p,
      const smt::SmtSolver & solver);
  ~Bmc();

  typedef Prover super;

  void initialize() override;

  ProverResult check_until(int k) override;

 private:
  bool step(int i);

};  // class Bmc

}  // namespace pono
