/*********************                                                        */
/*! \file term_analysis.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Useful functions for term analysis.
**
**
**/

#include "smt-switch/smt.h"

using namespace smt;

namespace pono {

void get_free_symbols(const Term & term, UnorderedTermSet & out_symbols)
{
  TermVec to_visit({ term });
  UnorderedTermSet visited;

  Term t;
  while (to_visit.size()) {
    t = to_visit.back();
    to_visit.pop_back();

    if (visited.find(t) == visited.end()) {
      visited.insert(t);
      // add children to queue
      for (auto tt : t) {
        to_visit.push_back(tt);
      }

      if (t->is_symbolic_const()) {
        out_symbols.insert(t);
      }
    }
  }
}

UnorderedTermSet get_free_symbols(const Term & term)
{
  UnorderedTermSet free_symbols;
  get_free_symbols(term, free_symbols);
  return free_symbols;
}
}  // namespace pono
