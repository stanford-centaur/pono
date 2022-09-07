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

#include "utils/logger.h"

using namespace smt;

namespace pono {

BmcSimplePath::BmcSimplePath(const Property & p, const TransitionSystem & ts,
                             const SmtSolver & solver,
                             PonoOptions opt)
  : super(p, ts, solver, opt)
{
  engine_ = Engine::BMC_SP;
  kind_engine_name_ = "BMC-SP";
  // BMC with simple path checking is a special instance of
  // k-induction; disable inductive case checking
  options_.kind_no_ind_check_ = true;
}

BmcSimplePath::~BmcSimplePath() {}

ProverResult BmcSimplePath::check_until(int k)
{
  initialize();
  super::check_until(k);
}

}  // namespace pono
