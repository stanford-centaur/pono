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
  FunctionalTransitionSystem(smt::SmtSolver & s) : TransitionSystem(s) {}

  /* deleted methods
     some methods are not supported for functional transition systems
     which ensure that the transition relation is composed entirely
     of invariants and equalities where the left side is a single
     next state symbol and the right is only over current states
     and inputs
  */
  void set_behavior(const smt::Term & init, const smt::Term & trans) override
  {
    throw CosaException("Can't call set_behavior on a FunctionalTransitionSystem");
  }
  void set_trans(const smt::Term & trans) override
  {
    throw CosaException("Can't call set_trans on a FunctionalTransitionSystem");
  }
  void constrain_trans(const smt::Term & constraint) override
  {
    throw CosaException("Can't call constrain_trans on a FunctionalTransitionSystem");
  }
  smt::Term next(const smt::Term & term) const override
  {
    throw CosaException("Can't call next on a FunctionalTransitionSystem");
  }
  bool is_next_var(const smt::Term & sv) const override
  {
    throw CosaException("Can't call set_behavior on a FunctionalTransitionSystem");
  }

  // overloaded
  smt::Term make_input(const std::string name, const smt::Sort & sort) override;

  // overloaded
  smt::Term make_state(const std::string name, const smt::Sort & sort) override;

  // overloaded
  bool is_functional() const override { return true; };

 protected:

  // helpers and checkers

  // overloaded
  bool known_symbols(const smt::Term & term);

  // Note: this method is not exposed, because the user should not be able to
  //       build terms with next states
  /* Replace all current states by their next-state updates, functionally */
  smt::Term to_next_func(const smt::Term & term);
};

}  // namespace cosa
