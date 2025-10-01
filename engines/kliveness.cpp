/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Po-Chun Chien
 ** This file is part of the pono project.
 ** Copyright (c) 2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief Implementation of k-liveness
 **
 **
 **/

#include "kliveness.h"

#include "utils/logger.h"

namespace pono {

KLiveness::KLiveness(const LivenessProperty & p,
                     const TransitionSystem & ts,
                     const smt::SmtSolver & solver,
                     PonoOptions opt)
    : super(p, ts, solver, opt)
{
  engine_ = opt.engine_;
  if (justice_conditions_.size() > 1) {
    throw PonoException("K-liveness only supports one justice condition.");
  }
}

KLiveness::~KLiveness() {}

void KLiveness::initialize()
{
  if (initialized_) {
    return;
  }
  super::initialize();
  solver_->reset_assertions();
}

ProverResult KLiveness::check_until(int k) { throw PonoException("NYI"); }

}  // namespace pono
