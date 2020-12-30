/*********************                                                        */
/*! \file make_provers.cpp
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

#include "make_provers.h"

#include "engines/bmc.h"
#include "engines/bmc_simplepath.h"
#include "engines/ic3ia.h"
#include "engines/interpolantmc.h"
#include "engines/kinduction.h"
#include "engines/mbic3.h"
#ifdef WITH_MSAT_IC3IA
#include "engines/msat_ic3ia.h"
#endif

#include "smt/available_solvers.h"

using namespace smt;
using namespace std;

namespace pono {

vector<Engine> all_engines() { return { BMC, BMC_SP, KIND, INTERP, MBIC3 }; }

shared_ptr<Prover> make_prover(Engine e,
                               Property & p,
                               SmtSolver & slv,
                               PonoOptions opts)
{
  if (e == BMC) {
    return make_shared<Bmc>(p, slv, opts);
  } else if (e == BMC_SP) {
    return make_shared<BmcSimplePath>(p, slv, opts);
  } else if (e == KIND) {
    return make_shared<KInduction>(p, slv, opts);
  } else if (e == INTERP) {
#ifdef WITH_MSAT
    SmtSolver s = create_interpolating_solver(SolverEnum::MSAT_INTERPOLATOR);
    return make_prover(e, p, slv, s, opts);
#else
    throw PonoException(
        "Interpolant-based modelchecking requires an interpolator");
#endif
  } else if (e == MBIC3) {
    return make_shared<ModelBasedIC3>(p, slv, opts);
  } else if (e == IC3IA_ENGINE) {
#ifdef WITH_MSAT
    SmtSolver s = create_interpolating_solver(SolverEnum::MSAT_INTERPOLATOR);
    return make_prover(e, p, slv, s, opts);
#else
    throw PonoException(
        "IC3IA uses MathSAT for interpolants, but not built with MathSAT");
#endif
#ifdef WITH_MSAT_IC3IA
  } else if (e == MSAT_IC3IA) {
    return make_shared<MsatIC3IA>(p, slv, opts);
#endif
  } else {
    throw PonoException("Unhandled engine");
  }
}

shared_ptr<Prover> make_prover(Engine e,
                               Property & p,
                               SmtSolver & slv,
                               SmtSolver & itp,
                               PonoOptions opts)
{
  if (e == INTERP) {
    return make_shared<InterpolantMC>(p, slv, itp, opts);
  } else if (e == IC3IA_ENGINE) {
    return make_shared<IC3IA>(p, slv, itp, opts);
  } else {
    throw PonoException(
        "Got unexpected engine when passing a solver and interpolator to "
        "make_prover.");
  }
}

}  // namespace pono
