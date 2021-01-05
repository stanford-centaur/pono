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

#include "options/options.h"
#include "smt-switch/smt.h"

namespace pono {

/** Creates an SmtSolver of the provided type
 *  @param se the SolverEnum to identify which type of solver
 *  @param logging whether or not to keep track of term DAG at smt-switch level
 *         defaults to false because generally slower
 *  @param set the incremental option for the solver
 *  @param set the procude-model option for the solver
 *  @return an SmtSolver
 */
smt::SmtSolver create_solver(smt::SolverEnum se, bool logging=false,
                             bool incremental=true, bool produce_model=true);

// same as create_solver but will set reasonable options
// for particular engines (mostly IC3-variants)
smt::SmtSolver create_solver_for(smt::SolverEnum se, Engine e, bool logging);

/** Creates an interpolating SmtSolver of the provided type */
smt::SmtSolver create_interpolating_solver(smt::SolverEnum se);

// same as create_interpolating_solver but will set reasonable options
// for particular engines (mostly IC3-variants)
smt::SmtSolver create_interpolating_solver_for(smt::SolverEnum se, Engine e);

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
