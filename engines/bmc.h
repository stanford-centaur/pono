/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Ahmed Irfan, Makai Mann, Florian Lonsing
 ** This file is part of the pono project.
 ** Copyright (c) 2019, 2021, 2022 by the authors listed in the file AUTHORS
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
  // BMC bound to start with (default: 0)
  unsigned int bound_start_;
  // Value by which to increase BMC bound (default: 1)
  unsigned int bound_step_;
  // Used in binary search for cex: number of times we called 'solver->push()'
  unsigned int bin_search_frames_;
  // Get an upper bound on the cex, which is located in interval '[lb,ub]'
  int bmc_interval_get_cex_ub(const int lb, const int ub);
  // Add negated bad state predicate for all bounds in interval '[start,end]'.
  // This way, we restrict the search space of the solver to disregard these
  // bounds when searching for a cex.
  void bmc_interval_block_cex_ub(const int start, const int end);
  // Run linear search for cex within interval '[reached_k_ + 1, upper_bound]'
  void find_shortest_cex_linear_search(const int upper_bound);
  // Run binary search for cex within interval '[reached_k_ + 1, upper_bound]'.
  bool find_shortest_cex_binary_search(const int upper_bound);
  // Run binary search for cex within interval '[reached_k_ + 1, upper_bound]'
  // with less incremental solver use.
  bool find_shortest_cex_binary_search_less_inc(const int upper_bound);
};  // class Bmc

}  // namespace pono
