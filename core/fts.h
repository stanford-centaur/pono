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

#include <unordered_map>

#include "ts.h"
#include "utils/exceptions.h"

namespace cosa {

class FunctionalTransitionSystem : public TransitionSystem
{
 public:
  FunctionalTransitionSystem(smt::SmtSolver s) : TransitionSystem(s) {}

  // overloaded
  bool is_functional() const override { return true; };

 protected:

  // helpers and checkers

  // overloaded
  bool known_symbols(const smt::Term & term) const override;

  // Note: this method is not exposed, because the user should not be able to
  //       build terms with next states
  /* Replace all current states by their next-state updates, functionally */
  smt::Term to_next_func(const smt::Term & term);
};

}  // namespace cosa
