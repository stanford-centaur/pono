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

#include "core/ts.h"

#include "smt-switch/smt.h"

namespace pono {

/*
 * unroller inspired by the ic3ia unroller:
 * https://es-static.fbk.eu/people/griggio/ic3ia/
 *
 */

class Unroller
{
 public:
  Unroller(const TransitionSystem & ts, const smt::SmtSolver & solver);
  ~Unroller();

  smt::Term at_time(const smt::Term & t, unsigned int k);

  smt::Term untime(const smt::Term & t) const;

  /** Returns the time of an unrolled variable
   *  example: get_var_time(x@4) = 4
   *  this only works for unrolled variables
   *  @param v the unrolled variable to get the time of
   *  @return a non-negative integer time
   */
  size_t get_var_time(const smt::Term & v) const;

 private:
  smt::Term var_at_time(const smt::Term & v, unsigned int k);
  smt::UnorderedTermMap & var_cache_at_time(unsigned int k);

  const TransitionSystem & ts_;
  const smt::SmtSolver solver_;

  typedef std::vector<smt::UnorderedTermMap> TimeCache;
  TimeCache time_cache_;
  TimeCache time_var_map_;
  smt::UnorderedTermMap untime_cache_;
  std::unordered_map<smt::Term, size_t> var_times_;

};  // class Unroller

}  // namespace pono
