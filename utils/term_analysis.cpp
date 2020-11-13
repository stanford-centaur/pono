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
#include "smt-switch/utils.h"
#include "utils/exceptions.h"

using namespace smt;
using namespace std;

namespace pono {

// set of boolean operators
// boolean terms with these operators are not predicates
unordered_set<PrimOp> boolops(
    { And,
      Or,
      Xor,
      Not,
      Implies,
      Iff,
      // Note: also including bit-vector operators for solvers that
      //       alias bool and bv of size 1
      //       should not make a difference for solvers that don't
      //       alias, because it will check if it's a boolean first.
      //       so this is just to make this method work for all solvers
      BVAnd,
      BVOr,
      BVXor,
      BVNand,
      BVNot });

// helper functions

/** Helper function for get_combinations. Builds the output vector recursively
 *  @param options a sequence of non-empty vectors containing each option for
 * this position
 *  @param output vector (should start empty) - gets incrementally populated
 * with a flattened sequence of vectors with each option
 *  @param i the index this recursive function is currently working on
 *  NOTE: this is not the most performant implementation, but in practice we
 * don't expect too many options so this is not expected to be a bottleneck
 */
void get_combinations_helper(const vector<TermVec> & options,
                             vector<TermVec> & out,
                             size_t i)
{
  if (i == options.size()) {
    return;
  }

  if (i == 0) {
    if (!out.empty()) {
      throw PonoException("Expecting empty vector");
    }
    // need to initialize the vector
    for (auto c : options[i]) {
      out.push_back({ c });
    }
  } else {
    size_t orig_out_size = out.size();
    for (size_t j = 0; j < orig_out_size; ++j) {
      if (options[i].empty()) {
        throw PonoException("Each element of options must be non-empty");
      }

      for (int k = options[i].size() - 1; k >= 0; --k) {
        // each element of out should be the same length
        // exactly i because that's how far into the options we've gotten
        assert(out[j].size() == i);
        if (k != 0) {
          // copy the vector and add it to out with this option appended
          TermVec outj = out[j];
          outj.push_back(options[i][k]);
          out.push_back(outj);
        } else {
          // for element 0, just add it to the existing vector instead of
          // creating a new one NOTE: this is the last option we're appending so
          // we don't need to keep the old
          //       version of this vector
          out[j].push_back(options[i][k]);
        }
      }
    }
  }

  // recursive call
  get_combinations_helper(options, out, i + 1);
}

/** Generate all combinations of the terms in the vector
 *  e.g. if the input is
 *  [[v, w],
 *    [x],
 *    [y, z]
 *  ]
 *  then the output should be [[v, x, y], [v, x, z], [w, x, y], [w, x, z]]
 *  @param options a sequence of non-empty vectors containing each option for
 * this position
 *  @param a flattened sequence of vectors with each option
 */
vector<TermVec> get_combinations(const vector<TermVec> & options)
{
  vector<TermVec> out;
  get_combinations_helper(options, out, 0);
  return out;
}

// end helper functions

UnorderedTermSet get_free_symbols(const Term & term)
{
  UnorderedTermSet free_symbols;
  get_free_symbols(term, free_symbols);
  return free_symbols;
}

void get_predicates(const SmtSolver & solver,
                    const Term & term,
                    UnorderedTermSet & out,
                    bool include_symbols)
{
  // NOTE: this is better than checking the SortKind of a sort
  //       some solvers alias sorts and might return a SortKind
  //       of BV (with width 1) for a boolean
  //       But, even those solvers will be consistent and a term
  //       t with boolean sort will satisfy
  //       t->get_sort() == boolsort
  //       even if
  //       t->get_sort()->get_sort_kind() != BOOL
  Sort boolsort = solver->make_sort(BOOL);

  TermVec to_visit({ term });
  UnorderedTermSet visited;

  Term t;
  while (to_visit.size()) {
    t = to_visit.back();
    assert(t);  // non-null term
    to_visit.pop_back();

    if (visited.find(t) == visited.end()) {
      visited.insert(t);

      TermVec children(t->begin(), t->end());
      // later we will want to know which children are ITEs (if any)
      unordered_set<size_t> ite_indices;
      // add children to stack
      Term c;
      for (size_t i = 0; i < children.size(); ++i) {
        c = children[i];
        to_visit.push_back(c);
        if (c->get_op() == Ite) {
          ite_indices.insert(i);
        }
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
      } else if (t->is_value()) {
        continue;
      }

      Op op = t->get_op();
      // no case in smt-switch (yet) where boolean term that is not
      // a symbolic const will have a null operator
      assert(!op.is_null());

      // special case for ITE children
      // Note: we're trying to never include an ITE in a predicate
      //       so if we get y = ite(x < 10, x+1, 0), we want to add
      //       y = x+1 and y = 0 as the predicates instead of the
      //       whole formula
      if (ite_indices.size()) {
        vector<TermVec> options;
        for (size_t i = 0; i < children.size(); ++i) {
          if (ite_indices.find(i) != ite_indices.end()) {
            TermVec ite_children(children[i]->begin(), children[i]->end());
            assert(ite_children.size() == 3);
            options.push_back({ ite_children[1], ite_children[2] });
          } else {
            options.push_back({ children[i] });
          }
        }
        assert(options.size() == children.size());
        // generate all combinations of options
        vector<TermVec> all_combinations = get_combinations(options);

        // then rebuild for each TermVec of children
        Term res;
        for (auto comb : all_combinations) {
          // construct a new term with the given combination of children
          assert(comb.size() == children.size());
          res = solver->make_term(op, comb);
          // add this term to the stack of terms to check for predicates
          to_visit.push_back(res);
        }
      } else if (boolops.find(op.prim_op) == boolops.end()) {
        // boolean terms that do not use a boolean combination operator are
        // predicates
        out.insert(t);
      }
    }
  }
}

}  // namespace pono
