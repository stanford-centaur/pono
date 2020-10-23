/*********************                                                  */
/*! \file ic3ia.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief IC3 via Implicit Predicate Abstraction (IC3IA) implementation
**        based on
**
**        IC3 Modulo Theories via Implicit Predicate Abstraction
**            -- Alessandro Cimatti, Alberto Griggio,
**               Sergio Mover, Stefano Tonetta
**
**        and the open source implementation:
**
**        https://es-static.fbk.eu/people/griggio/ic3ia/index.html
**/

#include "engines/ic3ia.h"

using namespace smt;
using namespace std;

namespace pono {

IC3IA::IC3IA(Property & p, SolverEnum se)
    : super(p, se), abs_ts_(ts_.solver()), ia_(ts_, abs_ts_)
{
  initialize();
}

IC3IA::IC3IA(Property & p, const SmtSolver & slv)
    : super(p, slv), abs_ts_(ts_.solver()), ia_(ts_, abs_ts_)
{
  initialize();
}

IC3IA::IC3IA(const PonoOptions & opt, Property & p, const SolverEnum se)
    : super(opt, p, se), abs_ts_(ts_.solver()), ia_(ts_, abs_ts_)
{
  initialize();
}

IC3IA::IC3IA(const PonoOptions & opt, Property & p, const SmtSolver & slv)
    : super(opt, p, slv), abs_ts_(ts_.solver()), ia_(ts_, abs_ts_)
{
  initialize();
}

// protected methods
void IC3IA::initialize()
{
  // TODO add the initial predicats
  throw PonoException("NYI");
}

void IC3IA::add_predicate(const Term & pred) { throw PonoException("NYI"); }

}  // namespace pono
