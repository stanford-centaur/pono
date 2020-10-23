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

#include "assert.h"
#include "smt-switch/smt.h"

using namespace smt;
using namespace std;

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
      // add children to stack
      for (auto tt : t) {
        to_visit.push_back(tt);
      }

      if (t->is_symbol()) {
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

void get_predicates(const Term & term,
                    const Sort & boolsort,
                    UnorderedTermSet & out,
                    bool include_symbols)
{
  TermVec to_visit({ term });
  UnorderedTermSet visited;

  // set of boolean operators
  // boolean terms with these operators are not predicates
  unordered_set<PrimOp> boolops({ And, Or, Xor, Not, Implies, Iff, Ite });

  Term t;
  while (to_visit.size()) {
    t = to_visit.back();
    assert(t);  // non-null term
    to_visit.pop_back();

    if (visited.find(t) == visited.end()) {
      visited.insert(t);
      // add children to stack
      for (auto tt : t) {
        to_visit.push_back(tt);
      }

      if (t->get_sort() != boolsort) {
        // not a candidate for predicates
        continue;
      }

      if (t->is_symbol()) {
        if (include_symbols) {
          out.insert(t);
        }
        continue;
      }

      Op op = t->get_op();
      // no case in smt-switch (yet) where boolean term that is not
      // a symbolic const will have a null operator
      assert(!op.is_null());

      if (boolops.find(op.prim_op) == boolops.end()) {
        // boolean terms that do not use a boolean combination operator are
        // predicates
        out.insert(t);
      }
    }
  }
}

}  // namespace pono
