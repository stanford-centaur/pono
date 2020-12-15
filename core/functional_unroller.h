/*********************                                                        */
/*! \file functional_unroller.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief An unroller implementation for functional transition systems that has
**        a configurable parameter for when to introduce new timed variables.
**
**
**/
#pragma once

#include "core/unroller.h"

namespace pono {
class FunctionalUnroller : public Unroller
{
 public:
  /** Instantiate a FunctionalUnroller
   *  @param ts the transition system to unroll
   *  @param solver the solver to use
   *  @param interval -- the interval to introduce new timed variables for state
   * vars input variables need to be unrolled at every step but state variables
   * can be substituted for directly in a functional system for: interval == 1 :
   * equivalent to an Unroller interval > 1  : introduces fresh timed variables
   * very interval steps interval == 0 : never introduces fresh timed variables
   * (pure functional substitution in unrolling)
   */
  FunctionalUnroller(const TransitionSystem & ts,
                     const smt::SmtSolver & solver,
                     size_t interval);
  ~FunctionalUnroller() {}

  typedef Unroller super;

 protected:
  size_t interval_;

  // overridden to use the interval_ parameter as described in constructor
  // documentation
  smt::UnorderedTermMap & var_cache_at_time(unsigned int k) override;
};
}  // namespace pono
