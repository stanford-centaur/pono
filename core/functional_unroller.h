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
   *         variables.
   *         Input variables need to be unrolled at every step but state
   *         variables can be substituted for directly in a functional system.
   *
   *         interval == 0 : never introduces fresh timed variables
   *         interval == 1 : equivalent to a regular Unroller
   *         interval > 1  : introduces fresh timed variables every interval
   *                           steps (pure functional substitution in
   *                           unrolling). E.g. if 4, then every 4 timesteps
   *                           it introduce new symbols.
   *
   *  for non-zero interval_ there are extra constraints that need to be added
   *  for an unrolling (that give fresh symbols a meaning)
   *  These are available through extra_constraints_at
   *
   */
  FunctionalUnroller(const TransitionSystem & ts,
                     size_t interval = 0,
                     const std::string & time_identifier = "@");

  ~FunctionalUnroller() {}

  typedef Unroller super;

  /** takes a term over current state and input variables
   *  and returns it at a given time with a functional unrolling
   *  NOTE: will throw exception if there's a next-state variable
   *  which makes no sense in a functional unrolling
   *  @param t the term to unroll (for FunctionalUnroller must
   *         only have current variables and inputs)
   *  @param k the time to unroll the term at
   *  @return the unrolled term
   */
  smt::Term at_time(const smt::Term & t, unsigned int k) override;

  /** Provides extra constraints for a functional unrolling
   *  with intermittent fresh symbols
   *  this is an attempt to deal with ITE explosion in deeply
   *  nested substitutions
   *  Thus, every interval_ steps a new symbol is introduced
   *  even for state variables.
   *  However, we then need to add equality constraints to
   *  give those symbols meaning
   */
  smt::Term extra_constraints_at(unsigned int k)
  {
    if (k >= extra_constraints_.size()) {
      throw PonoException("Haven't unrolled enough for extra_constraints at "
                          + std::to_string(k));
    }
    return extra_constraints_.at(k);
  }

 protected:
  size_t interval_;

  smt::TermVec extra_constraints_;

  // useful term
  smt::Term true_;

  // overridden to use the interval_ parameter as described in constructor
  // documentation
  smt::UnorderedTermMap & var_cache_at_time(unsigned int k) override;
};
}  // namespace pono
