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
  BTOR, CVC4,

#if WITH_MSAT
      MSAT,
#endif

#if WITH_YICES2
      YICES2,
#endif
});

SmtSolver create_solver(SolverEnum se, bool logging, bool incremental,
                        bool produce_model)
{
  SmtSolver s;
  switch (se) {
    case BTOR: {
      s = BoolectorSolverFactory::create(logging);
      break;
    }
    case CVC4: {
      s = CVC4SolverFactory::create(logging);
      break;
    }
#if WITH_MSAT
    case MSAT: {
      s = MsatSolverFactory::create(logging);
      break;
    }
#endif
#if WITH_YICES2
    case YICES2: {
      s = Yices2SolverFactory::create(logging);
      break;
    }
#endif
    default: {
      throw SmtException("Unhandled solver enum");
    }
  }

  s->set_opt("incremental", incremental ? "true" : "false");
  s->set_opt("produce-models", produce_model ? "true" : "false");

  return s;
}

SmtSolver create_interpolating_solver(SolverEnum se)
{
  switch (se) {
#if WITH_MSAT
    // for convenience -- accept any MSAT SolverEnum
    case MSAT:
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
