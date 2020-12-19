/*********************                                                  */
/*! \file ic3sa.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief IC3 with Syntax-Guided Abstraction based on
**
**        Model Checking of Verilog RTL Using IC3 with Syntax-Guided
*Abstraction.
**            -- Aman Goel, Karem Sakallah
**
**
**  within Pono, we are building on the bit-level IC3 instead of directly
**  on IC3Base, because a lot of the functionality is the same
**  In particular, we don't need to override inductive generalization
**
**/

#include "engines/ic3sa.h"

#include "assert.h"

using namespace smt;
using namespace std;

namespace pono {

IC3SA::IC3SA(Property & p, smt::SolverEnum se) : super(p, se) {}

IC3SA::IC3SA(Property & p, const smt::SmtSolver & s) : super(p, s) {}

IC3SA::IC3SA(const PonoOptions & opt, Property & p, smt::SolverEnum se)
    : super(opt, p, se)
{
}

IC3SA::IC3SA(const PonoOptions & opt, Property & p, const smt::SmtSolver & s)
    : super(opt, p, s)
{
}

IC3Formula IC3SA::get_model_ic3formula(TermVec * out_inputs,
                                       TermVec * out_nexts) const
{
  throw PonoException("IC3SA::get_model_ic3formula NYI");
}

bool IC3SA::ic3formula_check_valid(const IC3Formula & u) const
{
  throw PonoException("IC3SA::ic3formula_check_valid NYI");
}

IC3Formula IC3SA::generalize_predecessor(size_t i, const IC3Formula & c)
{
  throw PonoException("IC3SA::generalize_predecessor NYI");
}

void IC3SA::check_ts() const { throw PonoException("IC3SA::check_ts NYI"); }

}  // namespace pono
