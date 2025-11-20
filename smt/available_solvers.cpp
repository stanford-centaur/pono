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

#include <cassert>
#include <iostream>
#include <unordered_set>
#include <vector>

#include "smt-switch/logging_solver.h"
#include "smt-switch/printing_solver.h"
#include "utils/logger.h"

// these two always included
#include "smt-switch/bitwuzla_factory.h"
#include "smt-switch/cvc5_factory.h"

#ifdef WITH_BOOLECTOR
#include "smt-switch/boolector_factory.h"
#endif

#ifdef WITH_MSAT
#include "smt-switch/msat_factory.h"
// these are for setting specific options
// e.g. in create_solver_for
#include "smt-switch/msat_solver.h"
#include "smt/msat_options.h"
#endif

#ifdef WITH_YICES2
#include "smt-switch/yices2_factory.h"
#endif

#ifdef WITH_Z3
#include "smt-switch/z3_factory.h"
#endif

using namespace smt;

namespace pono {

// set of SMT-solver options that are not allowed to be set externally
const std::unordered_set<std::string> disallowed_smt_opts(
    { "incremental",
      "produce-interpolants",
      "produce-models",
      "produce-unsat-assumptions" });

inline bool is_allowed_smt_opt(const std::string & opt)
{
  return disallowed_smt_opts.find(opt) == disallowed_smt_opts.end();
}

inline void check_allowed_smt_opts(const StringMap & opts)
{
  for (const auto & optpair : opts) {
    if (!is_allowed_smt_opt(optpair.first)) {
      throw SmtException("SMT-solver option '" + optpair.first
                         + "' is not allowed to be set externally.");
    }
  }
}

// list of regular (non-interpolator) solver enums
const std::vector<SolverEnum> solver_enums({
    BZLA,
    CVC5,
#ifdef WITH_BOOLECTOR
    BTOR,
#endif
#ifdef WITH_MSAT
    MSAT,
#endif
#ifdef WITH_YICES2
    YICES2,
#endif
#ifdef WITH_Z3
    Z3,
#endif
});

// internal method for creating a particular solver
// doesn't set any options
SmtSolver create_solver_base(SolverEnum se, bool logging, bool printing = false)
{
  SmtSolver s;
  auto printing_style = DEFAULT_STYLE;
  switch (se) {
    case BZLA: {
      s = BitwuzlaSolverFactory::create(logging);
      break;
    }
#ifdef WITH_BOOLECTOR
    case BTOR: {
      s = BoolectorSolverFactory::create(logging);
      break;
    }
#endif
    case CVC5: {
      s = Cvc5SolverFactory::create(logging);
      printing_style = CVC5_STYLE;
      break;
    }
#ifdef WITH_MSAT
    case MSAT: {
      s = MsatSolverFactory::create(logging);
      printing_style = MSAT_STYLE;
      break;
    }
#endif
#ifdef WITH_YICES2
    case YICES2: {
      s = Yices2SolverFactory::create(logging);
      break;
    }
#endif
#ifdef WITH_Z3
    case Z3: {
      s = Z3SolverFactory::create(logging);
      break;
    }
#endif
    default: {
      throw SmtException("Unhandled solver enum");
    }
  }
  if (printing) {
    s = create_printing_solver(s, &std::cerr, printing_style);
  }

  return s;
}

SmtSolver create_solver(SolverEnum se,
                        bool logging,
                        bool incremental,
                        bool produce_model,
                        bool printing,
                        const StringMap & solver_opts)
{
  check_allowed_smt_opts(solver_opts);
  SmtSolver s = create_solver_base(se, logging, printing);

  s->set_opt("incremental", incremental ? "true" : "false");
  s->set_opt("produce-models", produce_model ? "true" : "false");
  if (se == BTOR) {
    // BTOR does not support reset-assertions by default without this.
    s->set_opt("base-context-1", "true");
  }
  for (const auto & optpair : solver_opts) {
    s->set_opt(optpair.first, optpair.second);
  }

  return s;
}

SmtSolver create_solver_for(SolverEnum se,
                            Engine e,
                            bool logging,
                            bool full_model,
                            bool printing,
                            const StringMap & solver_opts)
{
  check_allowed_smt_opts(solver_opts);
  SmtSolver s;
  bool ic3_engine = ic3_variants().find(e) != ic3_variants().end();
  if (e == IC3SA_ENGINE) {
    // IC3SA requires a full model
    full_model = true;
  }

  if (se == MSAT && ic3_engine) {
#ifdef WITH_MSAT
    // These will be managed by the solver object
    // don't need to destroy
    StringMap opts(solver_opts);  // user-specified options take precedence
    opts.emplace("model_generation", "true");
    if (!full_model && e == IC3IA_ENGINE) {
      // only need boolean model
      opts.emplace("bool_model_generation", "true");
      opts.emplace("model_generation", "false");
      // Reasoning from open-source IC3IA code (by Alberto Griggio):
      // Turn off propagation of toplevel information. This is just overhead in
      // an IC3 context (where the solver is called hundreds of thousands of
      // times). Moreover, using it makes "lightweight" model generation (see
      // below) not effective
      opts.emplace("preprocessor.toplevel_propagation", "false");
    }
    msat_config cfg = get_msat_config_for_ic3(opts);
    msat_env env = msat_create_env(cfg);
    s = std::make_shared<MsatSolver>(cfg, env);
    if (logging) {
      s = std::make_shared<LoggingSolver>(s);
    }
    if (printing) {
      s = create_printing_solver(s, &std::cerr, MSAT_STYLE);
    }
    return s;
#else
    throw PonoException("not built with MathSAT");
#endif
  } else {
    // no special options yet for solvers other than mathsat
    s = create_solver(se, logging, true, true, printing, solver_opts);
  }

  assert(s);
  if (ic3_engine) {
    s->set_opt("produce-unsat-assumptions", "true");
  }
  return s;
}

SmtSolver create_reducer_for(SolverEnum se,
                             Engine e,
                             bool logging,
                             const StringMap & solver_opts)
{
  check_allowed_smt_opts(solver_opts);
  SmtSolver s;
  if (se != MSAT) {
    s = create_solver_base(se, logging);
    s->set_opt("incremental", "true");
    s->set_opt("produce-unsat-assumptions", "true");
    for (const auto & optpair : solver_opts) {
      s->set_opt(optpair.first, optpair.second);
    }
  } else {
#ifdef WITH_MSAT
    // user-specified options take precedence
    StringMap opts(solver_opts);
    // no models needed for a reducer
    opts.emplace("model_generation", "false");
    msat_config cfg = get_msat_config_for_ic3(opts);
    msat_env env = msat_create_env(cfg);
    s = std::make_shared<MsatSolver>(cfg, env);
    if (logging) {
      s = std::make_shared<LoggingSolver>(s);
    }
#else
    throw PonoException("not built with MathSAT");
#endif
  }

  assert(s);
  return s;
}

SmtSolver create_interpolating_solver(SolverEnum se,
                                      bool printing,
                                      const StringMap & solver_opts)
{
  check_allowed_smt_opts(solver_opts);
  SmtSolver s;
  PrintingStyleEnum printing_style = DEFAULT_STYLE;
  switch (se) {
    case CVC5:
    case CVC5_INTERPOLATOR: {
      s = Cvc5SolverFactory::create_interpolating_solver();
      printing_style = CVC5_STYLE;
      break;
    }
#ifdef WITH_MSAT
    case MSAT:
    case MSAT_INTERPOLATOR: {
      s = MsatSolverFactory::create_interpolating_solver();
      printing_style = MSAT_STYLE;
      break;
    }
#endif
    default: {
      throw SmtException("Unhandled solver enum");
    }
  }
  if (printing) {
    s = create_printing_solver(s, &std::cerr, printing_style);
  }
  for (const auto & optpair : solver_opts) {
    s->set_opt(optpair.first, optpair.second);
  }
  return s;
}

SmtSolver create_interpolating_solver_for(SolverEnum se,
                                          Engine e,
                                          bool printing,
                                          const StringMap & solver_opts)
{
  if (ic3_variants().find(e) == ic3_variants().end() || se != MSAT) {
    return create_interpolating_solver(se, printing, solver_opts);
  }

  check_allowed_smt_opts(solver_opts);
  switch (se) {
#ifdef WITH_MSAT
      // for convenience -- accept any MSAT SolverEnum
    case MSAT:
    case MSAT_INTERPOLATOR: {
      // These will be managed by the solver object
      // don't need to destroy
      StringMap msat_opts(solver_opts);
      msat_opts.emplace("bool_model_generation", "false");
      msat_opts.emplace("model_generation", "true");
      msat_config cfg = get_msat_config_for_ic3(msat_opts);
      msat_env env = msat_create_env(cfg);
      SmtSolver solver = std::make_shared<MsatInterpolatingSolver>(cfg, env);
      if (printing) {
        solver = create_printing_solver(solver, &std::cerr, MSAT_STYLE);
      }
      return solver;
    }
#endif
    default: {
      throw SmtException("Unhandled solver enum");
    }
  }
}

SmtSolver create_quantifier_solver(SolverEnum se,
                                   bool logging,
                                   bool printing,
                                   const StringMap & solver_opts)
{
  check_allowed_smt_opts(solver_opts);
  const SolverEnum fallback_se = CVC5;
  const auto se_attribs = get_solver_attributes(se);
  if (se_attribs.find(QUANTIFIERS) == se_attribs.end()) {
    logger.log(1,
               "WARNING: Solver {} does not support quantifiers, "
               "using {} instead.",
               to_string(se),
               to_string(fallback_se));
    se = fallback_se;
  } else if (se == BTOR || se == MSAT) {
    // Boolector may run into segfaults; MathSAT oftentimes return UNKNOWN
    logger.log(
        1,
        "WARNING: Solver {} does not work well with quantifiers in practice, "
        "using {} instead.",
        to_string(se),
        to_string(fallback_se));
    se = fallback_se;
  }
  SmtSolver s = create_solver_base(se, logging, printing);
  for (const auto & optpair : solver_opts) {
    s->set_opt(optpair.first, optpair.second);
  }
  return s;
}

const std::vector<SolverEnum> itp_enums({
    CVC5_INTERPOLATOR,
#ifdef WITH_MSAT
    MSAT_INTERPOLATOR,
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
