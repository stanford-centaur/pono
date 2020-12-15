/*********************                                                        */
/*! \file functional_unroller.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief An unroller implementation for functional transition systems that has
**        a configurable parameter for when to introduce new timed variables.
**
**
**/
#include "core/unroller.h"
#include "utils/exceptions.h"

using namespace smt;
using namespace std;

namespace pono {

FunctionalUnroller::FunctionalUnroller(const TransitionSystem & ts,
                                       const SmtSolver & solver,
                                       size_t interval)
    : super(ts, solver), interval_(interval)
{
  if (!ts.is_functional()) {
    throw PonoException(
        "Can only use FunctionalUnroller on a FunctionalTransitionSystem.");
  }
}

UnorderedTermMap & FunctionalUnroller::var_cache_at_time(unsigned int k)
{
  throw PonoException("NYI");
}

}  // namespace pono
