/*********************                                                  */
/*! \file functionalize.h
** \verbatim
** Top contributors (to current version):
**   Po-Chun Chien
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Turns a relational transition system into a functional one
**
**
**/

#pragma once

#include "core/ts.h"

namespace pono {

/** Returns a functional transition system with the traces of a relational
 *  one.
 *
 *  Let T be the transition relation without the conjuncts that the
 *  constraints of the system added to it, as the result keeps them as its
 *  own constraints. The result tries to need as few auxiliary inputs as
 *  possible: a conjunct x' = f(X, I) of T, with the sides in either order,
 *  gives the state variable x the update x' := f instead of an input. Here f
 *  may also refer to the next-state variables of other state variables that
 *  get an update so, which are replaced by their updates.
 *
 *  The conjuncts R of T that define no update are encoded relationally, as a
 *  transition relation is encoded as a circuit (see Sect. 5.2.1 of
 *  https://doi.org/10.48550/arXiv.2608.07397).
 *  For each other state variable x whose next-state variable occurs in R, a
 *  fresh auxiliary input x_in chooses the next value, x' := x_in. A fresh state
 *  variable valid, which starts true, records whether every step so far kept
 *  to R:
 *    valid' := valid /\ R(X, I, X')
 *  where each next-state variable is replaced by its update. valid is added
 *  as an invariant constraint, so a trace of the result is a trace of the
 *  original system that may end in a state without a successor. Without R,
 *  there is no valid. A state variable whose next-state variable does not
 *  occur in T stays without an update.
 *
 *  The variables, initial states and constraints of the original system are
 *  kept, so a counterexample keeps its length and the values of the original
 *  variables. The new variables are named as pono's own.
 *
 *  @param ts the transition system to functionalize
 *  @param allow_new_inputvars whether each x_in is an input, or else a state
 *         variable without an update, which takes any value in each step as
 *         well, for engines that take no inputs (e.g., IC3SA)
 *  @return ts itself if it is functional already, the functionalized system
 *          otherwise
 */
TransitionSystem functionalize(const TransitionSystem & ts,
                               bool allow_new_inputvars = true);

}  // namespace pono
