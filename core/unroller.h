/*********************                                                        */
/*! \file 
 ** \verbatim
 ** Top contributors (to current version):
 **   Ahmed Irfan, Makai Mann
 ** This file is part of the cosa2 project.
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

#include "ts.h"

#include "smt-switch/smt.h"

namespace cosa {

/*
 * unroller inspired by the ic3ia unroller:
 * https://es-static.fbk.eu/people/griggio/ic3ia/
 *
 */

class Unroller
{
 public:
  Unroller(const TransitionSystem & ts, smt::SmtSolver & solver);
  ~Unroller();

  smt::Term at_time(const smt::Term & t, unsigned int k);

  smt::Term untime(const smt::Term & t) const;

 private:
  smt::Term var_at_time(const smt::Term & v, unsigned int k);
  smt::UnorderedTermMap & var_cache_at_time(unsigned int k);

  const TransitionSystem & ts_;
  smt::SmtSolver & solver_;

  typedef std::vector<smt::UnorderedTermMap> TimeCache;
  TimeCache time_cache_;
  TimeCache time_var_map_;
  smt::UnorderedTermMap untime_cache_;

};  // class Unroller

}  // namespace cosa
