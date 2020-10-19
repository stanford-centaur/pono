/*********************                                                        */
/*! \file msat_ic3ia.h
** \verbatim
** Top contributors (to current version):
**   Ahmed Irfan, Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief A backend using the open-source IC3IA implementation from
**        Alberto Griggio. Available here:
**        https://es-static.fbk.eu/people/griggio/ic3ia/index.html
**
**/

#pragma once

#include "engines/prover.h"

#include "ic3.h"

namespace pono {

class MsatIC3IA : public Prover
{
 public:
  MsatIC3IA(Property & p, smt::SolverEnum se);
  MsatIC3IA(Property & p, const smt::SmtSolver & solver);
  MsatIC3IA(const PonoOptions & opt, Property & p, smt::SolverEnum se);
  MsatIC3IA(const PonoOptions & opt, Property & p, const smt::SmtSolver & solver);

  typedef Prover super;

  ProverResult prove() override;

  // TODO: might have to ignore this and only support prove
  // ic3ia doesn't have a way to only check up to some number of frames
  // might need to make prove virtual as well
  ProverResult check_until(int k) override;

};

}  // namespace pono
