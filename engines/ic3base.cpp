/*********************                                                  */
/*! \file ic3base.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan, Florian Lonsing
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Abstract base class implementation of IC3 parameterized by
**        the unit used in frames, pre-image computation, and inductive
**        and predecessor generalization techniques.
**
**/

#include <type_traits>

#include "engines/ic3base.h"

using namespace smt;
using namespace std;

namespace pono {

IC3Base::IC3Base(Property & p, smt::SolverEnum se) : super(p, se)
{
  initialize();
}

IC3Base::IC3Base(Property & p, const smt::SmtSolver & s) : super(p, s)
{
  initialize();
}

IC3Base::IC3Base(const PonoOptions & opt, Property & p, smt::SolverEnum se)
    : super(opt, p, se)
{
  initialize();
}
IC3Base::IC3Base(const PonoOptions & opt,
                 Property & p,
                 const smt::SmtSolver & s)
    : super(opt, p, s)
{
  initialize();
}

IC3Base::initialize()
{
  static_assert(std::is_base_of<IC3Unit, Unit>::value,
                "Expecting Unit to be a derived class of IC3Unit");
}

template class IC3Base<IC3Unit>;

}  // namespace pono
