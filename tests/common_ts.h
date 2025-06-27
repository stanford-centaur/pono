/*********************                                                        */
/*! \file common_ts.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Getters for common transition systems used in tests
**
**
**/

#include "core/ts.h"
#include "smt-switch/smt.h"

namespace pono_tests {

/** Creates state variable x, initializes to zero and counts up until max_val
 *  and then resets
 *  @param ts the transition system to add to (assumed to be empty)
 *  @param max_val the value to count up to before resetting
 */
void counter_system(pono::TransitionSystem & ts, const smt::Term & max_val);

}  // namespace pono_tests
