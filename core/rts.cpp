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
  auto curr = curr_map_.find(term);
  if (curr != curr_map_.end()) {
    no_state_updates_.erase(curr->second);
  }
  for (const auto & subTerm : *term) {
    set_updated_states(subTerm);
  }
}
}  // namespace pono
