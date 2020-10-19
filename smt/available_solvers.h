/*********************                                                        */
/*! \file available_solvers.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Utility functions for creating solvers and checking for which
**        solvers this version of pono was built with.
**
**/

#pragma once

#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "smt-switch/smt.h"

// these two always included
#include "smt-switch/boolector_factory.h"
#include "smt-switch/cvc4_factory.h"

#if WITH_MSAT
#include "smt-switch/msat_factory.h"
#endif

#if WITH_YICES2
#include "smt-switch/yices2_factory.h"
#endif

namespace pono {

/** Creates an SmtSolver of the provided type
 *  @param se the SolverEnum to identify which type of solver
 *  @param logging whether or not to keep track of term DAG at smt-switch level
 *         defaults to false because generally slower
 *  @return an SmtSolver
 */
smt::SmtSolver create_solver(smt::SolverEnum se, bool logging=false);

/** Creates an interpolating SmtSolver of the provided type */
smt::SmtSolver create_interpolating_solver(smt::SolverEnum se);

// collect all the available solvers
std::vector<smt::SolverEnum> available_solver_enums();

// collect all the available interpolating solvers
std::vector<smt::SolverEnum> available_interpolator_enums();

/** Filter the available solvers by a set of attributes
 * @return all available solvers that have *all* the attributes
 */
std::vector<smt::SolverEnum> filter_solver_enums(
    const std::unordered_set<smt::SolverAttribute> attributes);

}  // namespace pono
