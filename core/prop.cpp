/*********************                                                        */
/*! \file 
 ** \verbatim
 ** Top contributors (to current version):
 **   Ahmed Irfan, Makai Mann
 ** This file is part of the cosa2 project.
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

namespace cosa {

Property::Property(const TransitionSystem & ts, Term p) : ts_(ts), prop_(p)
{
  const UnorderedTermSet & states = ts.states();
  UnorderedTermSet free_symbols = get_free_symbols(p);
  for (auto s : free_symbols) {
    if (states.find(s) == states.end()) {
      throw CosaException("Property should only use state variables");
    }
  }
}

Property::~Property() {}

}  // namespace cosa
