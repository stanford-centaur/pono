#pragma once

#include <iostream>
#include <unordered_map>
#include <vector>

#include "smt-switch/smt.h"

// Always include boolector
#include "smt-switch/boolector_factory.h"

#if WITH_CVC4
#include "smt-switch/cvc4_factory.h"
#endif

#if WITH_MSAT
#include "smt-switch/msat_factory.h"
#endif

namespace pono_tests {

typedef ::smt::SmtSolver (*create_solver_fun)(bool);

enum SolverEnum
{
  BTOR = 0,
  CVC4,
  MSAT
};

typedef std::unordered_map<SolverEnum, create_solver_fun> CreateSolverFunsMap;

// Create a map from enums to available solver creation functions
CreateSolverFunsMap available_solvers();

// collect all the available solvers
std::vector<SolverEnum> available_solver_enums();

std::ostream & operator<<(std::ostream & o, SolverEnum e);

}  // namespace pono_tests
