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
#include "engines/cegar_ops_uf.h"
#include "engines/cegar_values.h"
#include "engines/dual_approx_reach.h"
#include "engines/ic3bits.h"
#include "engines/ic3ia.h"
#include "engines/ic3sa.h"
#include "engines/interp_seq_mc.h"
#include "engines/interpolantmc.h"
#include "engines/kinduction.h"
#include "engines/mbic3.h"
#include "engines/syguspdr.h"
#ifdef WITH_MSAT_IC3IA
#include "engines/msat_ic3ia.h"
#endif
#include "smt-switch/smt.h"

using namespace smt;
using namespace std;

namespace pono {

vector<Engine> all_engines()
{
  return { BMC, BMC_SP, KIND, MBIC3, INTERP, ISMC, IC3IA_ENGINE, IC3SA_ENGINE };
}

shared_ptr<SafetyProver> make_prover(Engine e,
                                     const SafetyProperty & p,
                                     const TransitionSystem & ts,
                                     const SmtSolver & slv,
                                     PonoOptions opts)
{
  if (e == BMC) {
    return make_shared<Bmc>(p, ts, slv, opts);
  } else if (e == BMC_SP) {
    return make_shared<BmcSimplePath>(p, ts, slv, opts);
  } else if (e == DAR) {
    return make_shared<DualApproxReach>(p, ts, slv, opts);
  } else if (e == KIND) {
    return make_shared<KInduction>(p, ts, slv, opts);
  } else if (e == INTERP) {
    return make_shared<InterpolantMC>(p, ts, slv, opts);
  } else if (e == ISMC) {
    return make_shared<InterpSeqMC>(p, ts, slv, opts);
  } else if (e == MBIC3) {
    return make_shared<ModelBasedIC3>(p, ts, slv, opts);
  } else if (e == IC3_BITS) {
    return make_shared<IC3Bits>(p, ts, slv, opts);
  } else if (e == IC3IA_ENGINE) {
    return make_shared<IC3IA>(p, ts, slv, opts);
#ifdef WITH_MSAT_IC3IA
  } else if (e == MSAT_IC3IA) {
    return make_shared<MsatIC3IA>(p, ts, slv, opts);
#endif
  } else if (e == IC3SA_ENGINE) {
    return make_shared<IC3SA>(p, ts, slv, opts);
  } else if (e == SYGUS_PDR) {
    return make_shared<SygusPdr>(p, ts, slv, opts);
  } else {
    throw PonoException("Unhandled engine");
  }
}

shared_ptr<SafetyProver> make_ceg_proph_prover(Engine e,
                                               const SafetyProperty & p,
                                               const TransitionSystem & ts,
                                               const SmtSolver & slv,
                                               PonoOptions opts)
{
  if (e == BMC) {
    return std::make_shared<CegProphecyArrays<Bmc>>(p, ts, slv, opts);
  } else if (e == BMC_SP) {
    return std::make_shared<CegProphecyArrays<BmcSimplePath>>(p, ts, slv, opts);
  } else if (e == DAR) {
    return std::make_shared<CegProphecyArrays<DualApproxReach>>(
        p, ts, slv, opts);
  } else if (e == KIND) {
    return std::make_shared<CegProphecyArrays<KInduction>>(p, ts, slv, opts);
  } else if (e == INTERP) {
    return std::make_shared<CegProphecyArrays<InterpolantMC>>(p, ts, slv, opts);
  } else if (e == ISMC) {
    return std::make_shared<CegProphecyArrays<InterpSeqMC>>(p, ts, slv, opts);
  } else if (e == MBIC3) {
    return std::make_shared<CegProphecyArrays<ModelBasedIC3>>(p, ts, slv, opts);
  } else if (e == IC3IA_ENGINE) {
    return std::make_shared<CegProphecyArrays<IC3IA>>(p, ts, slv, opts);
  }
#ifdef WITH_MSAT_IC3IA
  else if (e == MSAT_IC3IA) {
    return std::make_shared<CegProphecyArrays<MsatIC3IA>>(p, ts, slv, opts);
  }
#endif
  else if (e == IC3SA_ENGINE) {
    return std::make_shared<CegProphecyArrays<IC3SA>>(p, ts, slv, opts);
  } else {
    throw PonoException("Unhandled engine");
  }
}

shared_ptr<SafetyProver> make_cegar_values_prover(Engine e,
                                                  const SafetyProperty & p,
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

template <class Prover_T>
shared_ptr<CegarOpsUf<Prover_T>> create_cegar_bv_arith_prover(
    const SafetyProperty & p,
    const TransitionSystem & ts,
    const SmtSolver & slv,
    PonoOptions opts,
    const UnorderedOpSet & ops_to_abstract)
{
  shared_ptr<CegarOpsUf<Prover_T>> prover =
      std::make_shared<CegarOpsUf<Prover_T>>(p, ts, slv, opts);
  prover->set_ops_to_abstract(ops_to_abstract);
  prover->set_min_bitwidth(opts.ceg_bv_arith_min_bw_);
  return prover;
}

shared_ptr<SafetyProver> make_cegar_bv_arith_prover(
    Engine e,
    const SafetyProperty & p,
    const TransitionSystem & ts,
    const SmtSolver & slv,
    PonoOptions opts,
    const UnorderedOpSet & ops_to_abstract)
{
  if (e == IC3IA_ENGINE) {
    if (opts.ceg_prophecy_arrays_) {
      return create_cegar_bv_arith_prover<CegProphecyArrays<IC3IA>>(
          p, ts, slv, opts, ops_to_abstract);
    } else {
      return create_cegar_bv_arith_prover<IC3IA>(
          p, ts, slv, opts, ops_to_abstract);
    }
  } else if (e == IC3SA_ENGINE) {
    return create_cegar_bv_arith_prover<IC3SA>(
        p, ts, slv, opts, ops_to_abstract);
  } else if (e == INTERP) {
    return create_cegar_bv_arith_prover<InterpolantMC>(
        p, ts, slv, opts, ops_to_abstract);
  } else if (e == KIND) {
    return create_cegar_bv_arith_prover<KInduction>(
        p, ts, slv, opts, ops_to_abstract);
  } else if (e == BMC) {
    return create_cegar_bv_arith_prover<Bmc>(p, ts, slv, opts, ops_to_abstract);
  } else if (e == BMC_SP) {
    return create_cegar_bv_arith_prover<BmcSimplePath>(
        p, ts, slv, opts, ops_to_abstract);
  } else {
    throw PonoException("CegarOpsUf does not support the chosen engine");
  }
}

}  // namespace pono
