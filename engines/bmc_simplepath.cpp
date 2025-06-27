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

#include "bmc_simplepath.h"

using namespace smt;

namespace pono {

BmcSimplePath::BmcSimplePath(const SafetyProperty & p,
                             const TransitionSystem & ts,
                             const SmtSolver & solver,
                             PonoOptions opt)
    : super(p, ts, solver, opt)
{
  engine_ = Engine::BMC_SP;
  kind_engine_name_ = "BMC-SP";

  // Options that affect k-induction also have an effect on BMC-SP

  // BMC with simple path checking is a special instance of
  // k-induction; disable inductive case checking based on property;
  // inductive case based on initial states is still checked as in
  // previous implementation of BMC-SP
  options_.kind_no_ind_check_property_ = true;
}

BmcSimplePath::~BmcSimplePath() {}

ProverResult BmcSimplePath::check_until(int k)
{
  initialize();
  return super::check_until(k);
}

}  // namespace pono
