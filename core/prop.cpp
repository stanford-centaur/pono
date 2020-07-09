/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Ahmed Irfan, Makai Mann
 ** This file is part of the pono project.
 ** Copyright (c) 2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief
 **
 **
 **/

#include "prop.h"
#include "exceptions.h"
#include "term_analysis.h"

using namespace smt;

namespace pono {

Property::Property(const TransitionSystem & ts, const Term & p)
    : ts_(ts), prop_(p)
{
  const UnorderedTermSet & states = ts.statevars();
  UnorderedTermSet free_symbols = get_free_symbols(p);
  for (auto s : free_symbols) {
    if (states.find(s) == states.end()) {
      throw PonoException("Property should only use state variables");
    }
  }
}

Property::Property(const Property & prop, TermTranslator & tt)
    : ts_(prop.ts_, tt),
      // only need to transfer if solvers are different
      // ts_ constructor does the same thing internally
      prop_((prop.transition_system().solver() == tt.get_solver())
                ? prop.prop_
                : tt.transfer_term(prop.prop_))
{
}

Property::~Property() {}

}  // namespace pono
