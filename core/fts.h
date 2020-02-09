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

#include "rts.h"

namespace cosa {

class FunctionalTransitionSystem : public RelationalTransitionSystem
{
 public:
 FunctionalTransitionSystem(smt::SmtSolver & s) : RelationalTransitionSystem(s)
  {
  }

  /* deleted methods
     some methods are not supported for functional transition systems
     which ensure that the transition relation is composed entirely
     of invariants and equalities where the left side is a single
     next state symbol and the right is only over current states
     and inputs
  */
  void set_behavior(const smt::Term & init, const smt::Term & trans) = delete;
  void set_trans(const smt::Term & trans) = delete;
  void constrain_trans(const smt::Term & constraint) = delete;
  smt::Term curr(const smt::Term & term) const = delete;
  smt::Term next(const smt::Term & term) const = delete;
  smt::Term is_curr_var(const smt::Term & sv) const = delete;
  smt::Term is_next_var(const smt::Term & sv) const = delete;

  // overloaded
  smt::Term make_input(const std::string name, const smt::Sort & sort);

  // overloaded
  smt::Term make_state(const std::string name, const smt::Sort & sort);


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
