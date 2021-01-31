/*********************                                                        */
/*! \file msat_options.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Special configuration options for mathsat
**
**
**/
#pragma once

#ifdef WITH_MSAT

#include <string>
#include <unordered_map>

#include "mathsat.h"
#include "utils/exceptions.h"

namespace pono {

using StringMap = std::unordered_map<std::string, std::string>;

// configuration options copied from open-source ic3ia implementation
// https://es-static.fbk.eu/people/griggio/ic3ia/index.html
msat_config get_msat_config_for_ic3(bool interp, const StringMap & options)
{
  msat_config cfg = msat_create_config();

  // no random decisions in the SAT solver
  msat_set_option(cfg, "dpll.branching_random_frequency", "0");

  // try not assigning values to theory atoms that occur only in
  // already-satisfied clauses. Example: given a clause (P | (t >= 0)), if P
  // is true, the value of (t >= 0) doesn't matter, and so it is a "ghost"
  msat_set_option(cfg, "dpll.ghost_filtering", "true");

  // Handling disequalities might be potentially quite expensive (especially
  // over the integers, where splitting of !(t = 0) into ((t < 0) | (t > 0))
  // is needed), so we want to avoid this as much as possible. If (t = 0)
  // occurs only positively in the input formula, but the SAT solver assigns
  // it to false, we can avoid sending the literal !(t = 0) to the
  // arithmetic solver, since if !(t = 0) causes an arithmetic conflict, we
  // can always flip it (as the atom never occurs negated in the input
  // formula, so assigning it to true can't cause any Boolean conflict)
  msat_set_option(cfg, "theory.la.pure_equality_filtering", "true");

  // Avoid splitting negated equalities !(t = 0) if t is of rational
  // type. Over the rationals, it is often cheaper to handle negated
  // equalities in a special way rather than relying on the general
  // splitting described above
  msat_set_option(cfg, "theory.la.split_rat_eq", "false");

  // Do not reset the internal scores for variables in the SAT solver
  // whenever a pop_backtrack_point() command is issued. In an IC3 context,
  // the solver is called very often on very similar problems, so it makes
  // sense to keep the variable scores computed so far, and maybe perform a
  // full solver reset only every few thousand calls rather than
  // reinitializing the scores at every pop()
  msat_set_option(cfg, "dpll.pop_btpoint_reset_var_order", "false");

  // enable interpolation if required
  msat_set_option(cfg, "interpolation", interp ? "true" : "false");
  msat_set_option(cfg, "preprocessor.interpolation_ite_elimination", "true");

  msat_set_option(cfg, "theory.bv.bit_blast_mode", "1");
  if (interp) {
    // interpolation for BV requires the lazy solver
    msat_set_option(cfg, "theory.bv.bit_blast_mode", "0");
    msat_set_option(cfg, "theory.bv.eager", "false");
  }

  for (const auto & optpair : options) {
    int fail =
        msat_set_option(cfg, optpair.first.c_str(), optpair.second.c_str());
    if (fail) {
      throw PonoException("Failed to set mathsat option: " + optpair.first
                          + " => " + optpair.second);
    }
  }

  return cfg;
}
}  // namespace pono

#endif
