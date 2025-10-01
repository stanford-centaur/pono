/*********************                                                        */
/*! \file rts.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Implements a relational transition system interface
**
**
**/

#include "rts.h"

#include "smt-switch/utils.h"

using namespace smt;
using namespace std;

namespace pono {

void RelationalTransitionSystem::set_behavior(const Term & init,
                                              const Term & trans)
{
  // TODO: Only do this check in debug mode
  if (!known_symbols(init) || !known_symbols(trans)) {
    throw PonoException("Unknown symbols");
  }
  init_ = init;
  trans_ = trans;
  set_updated_states(trans_);
}

void RelationalTransitionSystem::set_trans(const Term & trans)
{
  // TODO: Only do this check in debug mode
  if (!known_symbols(trans)) {
    throw PonoException("Unknown symbols");
  }
  trans_ = trans;
  set_updated_states(trans_);
}

void RelationalTransitionSystem::constrain_trans(const Term & constraint)
{
  // TODO: Only do this check in debug mode
  if (!known_symbols(constraint)) {
    throw PonoException("Unknown symbols");
  }
  trans_ = solver_->make_term(And, trans_, constraint);
  set_updated_states(constraint);
}

void RelationalTransitionSystem::set_updated_states(const smt::Term & term)
{
  UnorderedTermSet free_vars;
  get_free_symbolic_consts(term, free_vars);
  for (const auto & free_var : free_vars) {
    if (next_statevars_.find(free_var) != next_statevars_.end()) {
      no_state_updates_.erase(curr_map_[free_var]);
    }
  }
}
}  // namespace pono
