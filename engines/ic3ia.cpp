/*********************                                                  */
/*! \file ic3ia.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan
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
**
**  within Pono, we are building on the bit-level IC3 instead of directly
**  on IC3Base, because a lot of the functionality is the same
**  In particular, we don't need to override either of the generalization
**  functions. Instead focusing on abstract/refine.
**
**/

#include "engines/ic3ia.h"

namespace pono {

IC3IA::IC3IA(Property & p, smt::SolverEnum se) : super(p, se) {}

IC3IA::IC3IA(Property & p, const smt::SmtSolver & s) : super(p, s) {}

IC3IA::IC3IA(const PonoOptions & opt, Property & p, smt::SolverEnum se)
    : super(opt, p, se)
{
}

IC3IA::IC3IA(const PonoOptions & opt, Property & p, const smt::SmtSolver & s)
    : super(opt, p, s)
{
}

// pure virtual method implementations

IC3Formula IC3IA::get_ic3_formula() const { throw PonoException("NYI"); }

bool IC3IA::ic3_formula_check_valid(const IC3Formula & u) const
{
  throw PonoException("NYI");
}

void IC3IA::check_ts() const
{
  // basically a No-Op
  // no restrictions except that interpolants must be supported
  // instead of checking explicitly, just let the interpolator throw an
  // exception better than maintaining in two places
}

void IC3IA::initialize() { throw PonoException("NYI"); }

void IC3IA::abstract() { throw PonoException("NYI"); }

RefineResult IC3IA::refine() { throw PonoException("NYI"); }

}  // namespace pono
