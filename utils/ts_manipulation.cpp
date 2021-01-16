/*********************                                                        */
/*! \file ts_manipulation.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Useful functions for creating and manipulating transition systems
**
**
**/

#include "utils/ts_manipulation.h"

#include "core/fts.h"
#include "core/rts.h"

using namespace smt;
using namespace std;

namespace pono {

TransitionSystem create_fresh_ts(bool functional, const SmtSolver & solver)
{
  if (functional) {
    return FunctionalTransitionSystem(solver);
  } else {
    return RelationalTransitionSystem(solver);
  }
}

}  // namespace pono
