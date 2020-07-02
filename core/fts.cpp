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


#include "fts.h"

using namespace smt;
using namespace std;

namespace cosa {

// protected methods

bool FunctionalTransitionSystem::known_symbols(const Term & term) const
{
  return contains(term, UnorderedTermSetPtrVec{ &states_, &inputs_ });
}

Term FunctionalTransitionSystem::to_next_func(const Term & term)
{
  return solver_->substitute(term, state_updates_);
}

}  // namespace cosa
