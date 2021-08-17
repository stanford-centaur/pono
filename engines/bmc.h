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
  Bmc(const Property & p, const TransitionSystem & ts,
      const smt::SmtSolver & solver,
      PonoOptions opt = PonoOptions());

  ~Bmc();

  typedef Prover super;

  void initialize() override;

  ProverResult check_until(int k) override;

 protected:
  bool step(int i);

 private:
  int bmc_interval_get_cex_ub(const int lb, const int ub, const smt::Term bad_cl);
  void bmc_interval_find_shortest_cex(const int upper_bound);
  void bmc_interval_find_shortest_cex_binary_search(const int upper_bound); 
};  // class Bmc

}  // namespace pono
