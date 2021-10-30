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
#include "smt-switch/cvc5_factory.h"

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
  BTOR, CVC5,

#if WITH_MSAT
      MSAT,
#endif

#if WITH_YICES2
      YICES2,
#endif
});

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
    case CVC5: {
      s = Cvc5SolverFactory::create(logging);
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

SmtSolver create_solver_for(SolverEnum se,
                            Engine e,
                            bool logging,
                            bool full_model)
{
  SmtSolver s;
  bool ic3_engine = ic3_variants().find(e) != ic3_variants().end();
  if (e == IC3SA_ENGINE) {
    // IC3SA requires a full model
    full_model = true;
  }

  if (se != MSAT) {
    // no special options yet for solvers other than mathsat
    s = create_solver(se, logging);
  }
#ifdef WITH_MSAT
  else if (se == MSAT && ic3_engine) {
    // These will be managed by the solver object
    // don't need to destroy
    unordered_map<string, string> opts({ { "model_generation", "true" } });
    if (!full_model && e == IC3IA_ENGINE) {
      // only need boolean model
      opts["bool_model_generation"] = "true";
      opts["model_generation"] = "false";
      // Reasoning from open-source IC3IA code (by Alberto Griggio):
      // Turn off propagation of toplevel information. This is just overhead in
      // an IC3 context (where the solver is called hundreds of thousands of
      // times). Moreover, using it makes "lightweight" model generation (see
      // below) not effective
      opts["preprocessor.toplevel_propagation"] = "false";
    }
    msat_config cfg = get_msat_config_for_ic3(false, opts);
    msat_env env = msat_create_env(cfg);
    s = std::make_shared<MsatSolver>(cfg, env);
    if (logging) {
      s = make_shared<LoggingSolver>(s);
    }
    return s;
  }
#endif
  else {
    s = create_solver(se, logging);
  }

  assert(s);
  if (ic3_engine) {
    s->set_opt("produce-unsat-assumptions", "true");
  }
  return s;
}

SmtSolver create_reducer_for(SolverEnum se, Engine e, bool logging)
{
  SmtSolver s;
  if (se != MSAT) {
    s = create_solver_base(se, logging);
    s->set_opt("incremental", "true");
    s->set_opt("produce-unsat-assumptions", "true");
  }
#ifdef WITH_MSAT
  else {
    // no models needed for a reducer
    unordered_map<string, string> opts({ { "model_generation", "false" } });
    msat_config cfg = get_msat_config_for_ic3(false, opts);
    msat_env env = msat_create_env(cfg);
    s = std::make_shared<MsatSolver>(cfg, env);
    if (logging) {
      s = make_shared<LoggingSolver>(s);
    }
  }
#endif

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
  if (ic3_variants().find(e) == ic3_variants().end()) {
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
          get_msat_config_for_ic3(true,
                                  { { "bool_model_generation", "false" },
                                    { "model_generation", "true" } });
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

std::vector<SolverEnum> available_solver_enums_except(
    const std::unordered_set<SolverEnum> & exclude)
{
  std::vector<SolverEnum> res;
  for (const auto & se : solver_enums) {
    if (exclude.find(se) == exclude.end()) {
      res.push_back(se);
    }
  }
  return res;
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
