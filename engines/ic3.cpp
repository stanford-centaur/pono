/*********************                                                  */
/*! \file ic3.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Bit-level IC3 implementation using the IC3Base abstract base class
**/

#include "engines/ic3.h"

#include "assert.h"

using namespace smt;
using namespace std;

namespace pono {

// maps to expected operators for the two options for negated
// NOTE: uses both bv and boolean operators so that this works for Boolector
static const std::unordered_map<bool, std::unordered_set<smt::PrimOp>>
    expected_ops({ { false, { Or, BVOr } }, { true, { And, BVAnd } } });

// ClauseHandler implementation

IC3Unit ClauseHandler::create(const smt::TermVec & c) const
{
  assert(c.size());
  Term term = c.at(0);
  for (size_t i = 1; i < c.size(); ++i) {
    term = solver_->make_term(Or, term, c[i]);
  }
  IC3Unit res(term, c, false);
  assert(check_valid(res));
  return res;
}

IC3Unit ClauseHandler::negate(const IC3Unit & u) const
{
  const TermVec & children = u.children;
  assert(!u.is_null());
  assert(children.size());

  TermVec neg_children;
  neg_children.reserve(children.size());
  Term nc = smart_not(children.at(0));

  bool is_cube = u.negated;
  Term term = nc;
  neg_children.push_back(nc);
  for (size_t i = 1; i < children.size(); ++i) {
    nc = smart_not(children[i]);
    neg_children.push_back(nc);
    if (is_cube) {
      // negation is a clause
      term = solver_->make_term(Or, term, nc);
    } else {
      // negation is a cube
      term = solver_->make_term(And, term, nc);
    }
  }
  IC3Unit res(term, neg_children, !is_cube);
  return res;
}

bool ClauseHandler::check_valid(const IC3Unit & u) const
{
  bool is_valid = true;

  // check that children are literals
  Op op;
  for (auto c : u.children) {
    op = c->get_op();
    if (op == Not) {
      is_valid &= (*(c->begin()))->is_symbolic_const();
    } else {
      is_valid &= c->is_symbolic_const();
    }

    if (!is_valid) {
      return false;
    }
  }

  const unordered_set<PrimOp> & ops = expected_ops.at(u.negated);
  TermVec to_visit({ u.term });
  while (to_visit.size()) {
    Term t = to_visit.back();
    to_visit.pop_back();

    op = t->get_op();
    assert(!op.is_null());
    if (ops.find(op.prim_op) == ops.end()) {
      return false;
    }

    for (auto tt : t) {
      if (t->is_symbolic_const()
          || (t->get_op() == Not && (*(t->begin()))->is_symbolic_const())) {
        // hit a literal
        continue;
      }
      to_visit.push_back(tt);
    }
  }

  // got through all checks without failing
  assert(is_valid);
  return true;
}

}  // namespace pono
