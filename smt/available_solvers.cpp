/*********************                                                        */
/*! \file available_solvers.cpp
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

#include "available_solvers.h"

using namespace smt;

namespace pono {

// list of regular (non-interpolator) solver enums
const std::vector<SolverEnum> solver_enums({
  BTOR, BTOR_LOGGING, CVC4, CVC4_LOGGING,

#if WITH_MSAT
      MSAT, MSAT_LOGGING,
#endif

#if WITH_YICES2
      YICES2, YICES2_LOGGING,
#endif
});

SmtSolver create_solver(SolverEnum se)
{
  switch (se) {
    case BTOR: {
      return BoolectorSolverFactory::create(false);
      break;
      ;
    }
    case BTOR_LOGGING: {
      return BoolectorSolverFactory::create(true);
      break;
      ;
    }
    case CVC4: {
      return CVC4SolverFactory::create(false);
      break;
      ;
    }
    case CVC4_LOGGING: {
      return CVC4SolverFactory::create(true);
      break;
      ;
    }
#if WITH_MSAT
    case MSAT: {
      return MsatSolverFactory::create(false);
      break;
      ;
    }
    case MSAT_LOGGING: {
      return MsatSolverFactory::create(true);
      break;
      ;
    }
#endif
#if WITH_YICES2
    case YICES2: {
      return Yices2SolverFactory::create(false);
      break;
      ;
    }
    case YICES2_LOGGING: {
      return Yices2SolverFactory::create(true);
      break;
      ;
    }
#endif
    default: {
      throw SmtException("Unhandled solver enum");
    }
  }
}

SmtSolver create_interpolating_solver(SolverEnum se)
{
  switch (se) {
#if WITH_MSAT
    case MSAT_INTERPOLATOR: {
      return MsatSolverFactory::create_interpolating_solver();
      break;
      ;
    }
#endif
    default: {
      throw SmtException("Unhandled solver enum");
    }
  }
}

const std::vector<SolverEnum> itp_enums({
#if WITH_MSAT
  MSAT_INTERPOLATOR
#endif
});

std::vector<SolverEnum> available_solver_enums() { return solver_enums; }

std::vector<SolverEnum> available_no_logging_solver_enums()
{
  std::vector<SolverEnum> enums;
  for (auto se : solver_enums) {
    if (!is_logging_solver_enum(se)) {
      enums.push_back(se);
    }
  }
  return enums;
}

std::vector<SolverEnum> available_interpolator_enums() { return itp_enums; };

std::vector<SolverEnum> filter_solver_enums(
    const std::unordered_set<SolverAttribute> attributes)
{
  std::vector<SolverEnum> filtered_enums;
  for (auto se : solver_enums) {
    const std::unordered_set<SolverAttribute> & se_attrs =
        get_solver_attributes(se);

    bool all_attrs = true;
    for (auto a : attributes) {
      if (se_attrs.find(a) == se_attrs.end()) {
        all_attrs = false;
        break;
      }
    }

    if (all_attrs) {
      filtered_enums.push_back(se);
    }
  }

  return filtered_enums;
}

}  // namespace pono
