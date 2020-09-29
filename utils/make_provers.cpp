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
#include "engines/interpolantmc.h"
#include "engines/kinduction.h"
#include "engines/mbic3.h"

using namespace smt;
using namespace std;

namespace pono {

vector<Engine> all_engines() { return { BMC, BMC_SP, KIND, INTERP }; }

shared_ptr<Prover> make_prover(Engine e,
                               Property & p,
                               SolverEnum se,
                               PonoOptions opts)
{
  if (e == BMC) {
    return make_shared<Bmc>(opts, p, se);
  } else if (e == BMC_SP) {
    return make_shared<BmcSimplePath>(opts, p, se);
  } else if (e == KIND) {
    return make_shared<KInduction>(opts, p, se);
  } else if (e == INTERP) {
    return make_shared<InterpolantMC>(opts, p, se);
  } else if (e == MBIC3) {
    return make_shared<ModelBasedIC3>(opts, p, se);
  } else {
    throw PonoException("Unhandled engine");
  }
}

shared_ptr<Prover> make_prover(Engine e,
                               Property & p,
                               SmtSolver & slv,
                               PonoOptions opts)
{
  if (e == BMC) {
    return make_shared<Bmc>(opts, p, slv);
  } else if (e == BMC_SP) {
    return make_shared<BmcSimplePath>(opts, p, slv);
  } else if (e == KIND) {
    return make_shared<KInduction>(opts, p, slv);
  } else if (e == INTERP) {
    throw PonoException(
        "Interpolant-based modelchecking requires an interpolator");
  } else if (e == MBIC3) {
    return make_shared<ModelBasedIC3>(opts, p, slv);
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
    return make_shared<InterpolantMC>(opts, p, slv, itp);
  } else {
    throw PonoException(
        "Got unexpected engine when passing a solver and interpolator to "
        "make_prover.");
  }
}

}  // namespace pono
