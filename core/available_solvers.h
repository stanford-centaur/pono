#pragma once

#include <iostream>
#include <unordered_map>
#include <vector>

#include "smt-switch/smt.h"

namespace cosa {

/** Creates an SmtSolver of the provided type */
smt::SmtSolver create_solver(smt::SolverEnum se);

/** Creates an interpolating SmtSolver of the provided type */
smt::SmtSolver create_interpolator(smt::SolverEnum se);

// collect all the available solvers
std::vector<smt::SolverEnum> available_solver_enums();

// collect all the available interpolators
std::vector<smt::SolverEnum> available_interpolator_enums();

}  // namespace cosa
