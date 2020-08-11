/*********************                                                  */
/*! \file history_modifier.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Class for adding history variables that remember a value
**        for a certain number of steps to a system.
**
**/
#pragma once

#include "core/ts.h"

namespace pono {

class HistoryModifier
{
 public:
  HistoryModifier(TransitionSystem & ts);

  /** Returns a history variable with the given delay
   *  will create new variables and update the transition system as needed
   *  @param target a current state variable to target for a history variable
   *  @param delay the amount of delay to introduce for the history
   *  @return the history variable
   */
  smt::Term get_hist(const smt::Term & target, size_t delay);

 protected:
  TransitionSystem & ts_;
  const smt::SmtSolver solver_;

  // maps current state variables to a list of history variables
  // where the index corresponds to the delay (-1)
  std::unordered_map<smt::Term, smt::TermVec> hist_vars_;
};

}  // namespace pono
