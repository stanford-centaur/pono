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

Property::Property(const TransitionSystem & ts, const Term & p, std::string name)
  : ts_(ts), prop_(p), name_(name)
{
  const UnorderedTermSet & states = ts.statevars();
  UnorderedTermSet free_symbols = get_free_symbols(p);
  for (auto s : free_symbols) {
    if (!ts.is_curr_var(s))
    {
      throw PonoException("Property should only use state variables");
    }
  }

  // find a name if it wasn't provided
  // if no name is associated with it in the ts, then it will just
  // be the to_string of the term
  if (name_.empty())
  {
    name_ = ts_.get_name(prop_);
  }
}

Property::Property(const Property & prop, TermTranslator & tt)
    : ts_(prop.ts_, tt),
      // only need to transfer if solvers are different
      // ts_ constructor does the same thing internally
      prop_((prop.transition_system().solver() == tt.get_solver())
                ? prop.prop_
                : tt.transfer_term(prop.prop_)),
      name_(prop.name_)
{
}

Property::~Property() {}

}  // namespace pono
