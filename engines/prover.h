/*********************                                                        */
/*! \file 
 ** \verbatim
 ** Top contributors (to current version):
 **   Makai Mann, Ahmed Irfan
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

#include "prop.h"
#include "proverresult.h"
#include "smt-switch/smt.h"
#include "ts.h"
#include "unroller.h"

namespace cosa {
class Prover
{
 public:
  Prover(const Property & p, smt::SmtSolver & s);
  virtual ~Prover();

  virtual void initialize();

  virtual ProverResult check_until(int k) = 0;

  bool witness(std::vector<smt::UnorderedTermMap> & out);

  ProverResult prove();

 protected:
  const TransitionSystem & ts_;
  const Property & property_;

  smt::SmtSolver & solver_;
  Unroller unroller_;

  int reached_k_;

  smt::Term bad_;
};
}  // namespace cosa
