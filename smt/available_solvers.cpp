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

#include "smt/available_solvers.h"

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "assert.h"

// these two always included
#include "smt-switch/boolector_factory.h"
#include "smt-switch/cvc4_factory.h"

#if WITH_MSAT
#include "smt-switch/msat_factory.h"
// these are for setting specific options
// e.g. in create_solver_for
#include "smt-switch/logging_solver.h"
#include "smt-switch/msat_solver.h"
#include "smt/msat_options.h"

#endif

#if WITH_YICES2
#include "smt-switch/yices2_factory.h"
#endif

using namespace smt;
using namespace std;

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

// keep this up-to-date for setting solver options
// IC3 uses the solver in a different way, so different
// options are appropriate than for other engines
std::unordered_set<Engine> ic3_variants({ MBIC3, IC3IA_ENGINE, MSAT_IC3IA });

// internal method for creating a particular solver
// doesn't set any options
SmtSolver create_solver_base(SolverEnum se, bool logging)
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

  return s;
}

SmtSolver create_solver(SolverEnum se,
                        bool logging,
                        bool incremental,
                        bool produce_model)
{
  SmtSolver s = create_solver_base(se, logging);

  s->set_opt("incremental", incremental ? "true" : "false");
  s->set_opt("produce-models", produce_model ? "true" : "false");

  return s;
}

SmtSolver create_solver_for(SolverEnum se, Engine e, bool logging)
{
  if (se != MSAT && se != BTOR) {
    // no special options yet for solvers other than mathsat
    return create_solver(se, logging);
  }

  bool ic3_engine = ic3_variants.find(e) != ic3_variants.end();
  // special cases
  if (se == BTOR && e == IC3IA_ENGINE) {
    // for IC3IA it's best to be able to reset the solver
    // and boolector will do substitutions when there
    // are assertions at the base level
    // e.g. pred1 <-> p(X, Y)
    // then pred1 will be substituted for and no longer be
    // a symbol which causes problems for substitution, etc.
    // TODO adjust this based on whether we settle on using
    // variables for predicates in ic3ia
    SmtSolver s = create_solver(se, logging);
    s->set_opt("base-context-1", "true");
    return s;
  }
#ifdef WITH_MSAT
  else if (se == MSAT && ic3_engine) {
    // These will be managed by the solver object
    // don't need to destroy
    unordered_map<string, string> opts({ { "model_generation", "true" } });
    if (e == IC3IA_ENGINE) {
      // only need boolean model
      opts["bool_model_generation"] = "true";
      opts["model_generation"] = "false";
    }
    msat_config cfg = get_msat_config_for_ic3(false, opts);
    msat_env env = msat_create_env(cfg);
    SmtSolver s = std::make_shared<MsatSolver>(cfg, env);
    if (logging) {
      s = make_shared<LoggingSolver>(s);
    }
    return s;
  }
#endif
  else {
    return create_solver(se, logging);
  }
}

SmtSolver create_reducer_for(SolverEnum se, Engine e, bool logging)
{
  SmtSolver s;
  if (se == MSAT) {
    // no models needed for a reducer
    unordered_map<string, string> opts({ { "model_generation", "false" } });
    msat_config cfg = get_msat_config_for_ic3(false, opts);
    msat_env env = msat_create_env(cfg);
    s = std::make_shared<MsatSolver>(cfg, env);
    if (logging) {
      s = make_shared<LoggingSolver>(s);
    }
  } else {
    s = create_solver_base(se, logging);
    s->set_opt("incremental", "true");
    s->set_opt("produce-unsat-cores", "true");
  }

  if (se == BTOR && e == IC3IA_ENGINE) {
    s->set_opt("base-context-1", "true");
  }

  assert(s);
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

SmtSolver create_interpolating_solver_for(SolverEnum se, Engine e)
{
  if (ic3_variants.find(e) == ic3_variants.end()) {
    return create_interpolating_solver(se);
  }

  switch (se) {
#if WITH_MSAT
      // for convenience -- accept any MSAT SolverEnum
    case MSAT:
    case MSAT_INTERPOLATOR: {
      // These will be managed by the solver object
      // don't need to destroy
      msat_config cfg =
          get_msat_config_for_ic3(true, { { "model_generation", "false" } });
      msat_env env = msat_create_env(cfg);
      return std::make_shared<MsatInterpolatingSolver>(cfg, env);
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
