/*********************                                                        */
/*! \file make_provers.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Utility functions for creating provers.
**
**
**/

#pragma once

#include "engines/bmc.h"
#include "engines/bmc_simplepath.h"
#include "engines/interpolantmc.h"
#include "engines/kinduction.h"

namespace pono {

std::vector<Engine> all_engines() { return { BMC, BMC_SP, KIND, INTERP }; }

std::shared_ptr<Prover> make_prover(Engine e,
                                    const Property & p,
                                    smt::SolverEnum se,
                                    PonoOptions opts = PonoOptions())
{
  if (e == BMC) {
    return std::make_shared<Bmc>(opts, p, se);
  } else if (e == BMC_SP) {
    return std::make_shared<BmcSimplePath>(opts, p, se);
  } else if (e == KIND) {
    return std::make_shared<KInduction>(opts, p, se);
  } else if (e == INTERP) {
    return std::make_shared<InterpolantMC>(opts, p, se);
  } else {
    throw PonoException("Unhandled engine");
  }
}

std::shared_ptr<Prover> make_prover(Engine e,
                                    const Property & p,
                                    smt::SmtSolver & slv,
                                    PonoOptions opts = PonoOptions())
{
  if (e == BMC) {
    return std::make_shared<Bmc>(opts, p, slv);
  } else if (e == BMC_SP) {
    return std::make_shared<BmcSimplePath>(opts, p, slv);
  } else if (e == KIND) {
    return std::make_shared<KInduction>(opts, p, slv);
  } else if (e == INTERP) {
    throw PonoException(
        "Interpolant-based modelchecking requires an interpolator");
  } else {
    throw PonoException("Unhandled engine");
  }
}

std::shared_ptr<Prover> make_prover(Engine e,
                                    const Property & p,
                                    smt::SmtSolver & slv,
                                    smt::SmtSolver & itp,
                                    PonoOptions opts = PonoOptions())
{
  if (e == INTERP) {
    return std::make_shared<InterpolantMC>(opts, p, slv, itp);
  } else {
    throw PonoException(
        "Got unexpected engine when passing a solver and interpolator to "
        "make_prover.");
  }
}

}  // namespace pono
