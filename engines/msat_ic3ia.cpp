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

#include "engines/msat_ic3ia.h"

using namespace ic3ia;
using namespace smt;
using namespace std;

namespace pono {

MsatIC3IA::MsatIC3IA(Property & p, smt::SolverEnum se) : super(p, se)
{
  initialize();
}

MsatIC3IA::MsatIC3IA(Property & p, const SmtSolver & solver) : super(p, solver)
{
  initialize();
}

void MsatIC3IA::initialize() { throw PonoException("NYI"); }

ProverResult MsatIC3IA::check_until(int k) { throw PonoException("NYI"); }

}  // namespace pono
