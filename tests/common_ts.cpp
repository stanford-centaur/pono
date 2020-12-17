/*********************                                                        */
/*! \file common_ts.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Getters for common transition systems used in tests
**
**
**/

#include "tests/common_ts.h"

#include "assert.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

void counter_system(TransitionSystem & ts, const Term & max_val)
{
  assert(max_val);
  Sort sort = max_val->get_sort();
  Term x = ts.make_statevar("x", sort);
  SortKind sk = sort->get_sort_kind();
  PrimOp plus_op = (sk == BV) ? BVAdd : Plus;
  PrimOp lt_op = (sk == BV) ? BVUlt : Lt;
  Term inc_term = ts.make_term(plus_op, x, ts.make_term(1, sort));
  Term zero = ts.make_term(0, sort);
  ts.assign_next(
      x, ts.make_term(Ite, ts.make_term(lt_op, x, max_val), inc_term, zero));
  ts.set_init(ts.make_term(Equal, x, zero));
}

}  // namespace pono_tests
