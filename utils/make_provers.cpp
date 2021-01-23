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
#include "engines/ceg_prophecy_arrays.h"
#include "engines/cegar_values.h"
#include "engines/ic3bits.h"
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
                               const Property & p,
                               const TransitionSystem & ts,
                               const SmtSolver & slv,
                               PonoOptions opts)
{
  if (e == BMC) {
    return make_shared<Bmc>(p, ts, slv, opts);
  } else if (e == BMC_SP) {
    return make_shared<BmcSimplePath>(p, ts, slv, opts);
  } else if (e == KIND) {
    return make_shared<KInduction>(p, ts, slv, opts);
  } else if (e == INTERP) {
#ifdef WITH_MSAT
    return make_shared<InterpolantMC>(p, ts, slv, opts);
#else
    throw PonoException(
        "Interpolant-based modelchecking requires an interpolator");
#endif
  } else if (e == MBIC3) {
    return make_shared<ModelBasedIC3>(p, ts, slv, opts);
  } else if (e == IC3_BITS) {
    return make_shared<IC3Bits>(p, ts, slv, opts);
  } else if (e == IC3IA_ENGINE) {
#ifdef WITH_MSAT
    return make_shared<IC3IA>(p, ts, slv, opts);
#else
    throw PonoException(
        "IC3IA uses MathSAT for interpolants, but not built with MathSAT");
#endif
#ifdef WITH_MSAT_IC3IA
  } else if (e == MSAT_IC3IA) {
    return make_shared<MsatIC3IA>(p, ts, slv, opts);
#endif
  } else {
    throw PonoException("Unhandled engine");
  }
}

shared_ptr<Prover> make_ceg_proph_prover(Engine e,
                                         const Property & p,
                                         const TransitionSystem & ts,
                                         const SmtSolver & slv,
                                         PonoOptions opts)
{
  if (e == BMC) {
    return std::make_shared<CegProphecyArrays<Bmc>>(p, ts, slv, opts);
  } else if (e == BMC_SP) {
    return std::make_shared<CegProphecyArrays<BmcSimplePath>>(p, ts, slv, opts);
  } else if (e == KIND) {
    return std::make_shared<CegProphecyArrays<KInduction>>(p, ts, slv, opts);
  } else if (e == INTERP) {
    return std::make_shared<CegProphecyArrays<InterpolantMC>>(p, ts, slv, opts);
  } else if (e == MBIC3) {
    return std::make_shared<CegProphecyArrays<ModelBasedIC3>>(p, ts, slv, opts);
  } else if (e == IC3IA_ENGINE) {
#ifdef WITH_MSAT
    return std::make_shared<CegProphecyArrays<IC3IA>>(p, ts, slv, opts);
#else
    throw PonoException(
        "IC3IA uses MathSAT for interpolants, but not built with MathSAT");
#endif
  }
#ifdef WITH_MSAT_IC3IA
  else if (e == MSAT_IC3IA) {
    return std::make_shared<CegProphecyArrays<MsatIC3IA>>(p, ts, slv, opts);
  }
#endif
  else {
    throw PonoException("Unhandled engine");
  }
}

shared_ptr<Prover> make_cegar_values_prover(Engine e,
                                            const Property & p,
                                            const TransitionSystem & ts,
                                            const SmtSolver & slv,
                                            PonoOptions opts)
{
  if (e != IC3IA_ENGINE || !opts.ceg_prophecy_arrays_) {
    throw PonoException(
        "CegarValues currently only supports IC3IA with CegProphecyArrays");
  }

  return make_shared<CegarValues<CegProphecyArrays<IC3IA>>>(p, ts, slv, opts);
}

}  // namespace pono
