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

#include "rts.h"
#include "smt-switch/smt.h"

namespace cosa {

class Property
{
 public:
  Property(const TransitionSystem & ts, smt::Term p);
  ~Property();

  const smt::Term prop() const { return prop_; }

  const TransitionSystem & transition_system() const { return ts_; }

 private:
  const TransitionSystem & ts_;

  smt::Term prop_;

};  // class Property

}  // namespace cosa
