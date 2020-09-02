/*********************                                                        */
/*! \file adaptive_unroller.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief An unroller implementation that will work even when new variables
**        are added to the transition system on the fly.
**
**/

#pragma once

#include "core/unroller.h"

namespace pono {

class AdaptiveUnroller : public Unroller
{
 public:
  AdaptiveUnroller(const TransitionSystem & ts, const smt::SmtSolver & solver);

  typedef Unroller super;

 protected:
  smt::UnorderedTermMap & var_cache_at_time(unsigned int k) override;

  size_t num_vars_;  ///< the last known number of variables in the transition
                     ///< system
};

}  // namespace pono
