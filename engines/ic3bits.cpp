/*********************                                                  */
/*! \file ic3bits.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Bit-level IC3 implementation that splits bitvector variables
**        into the individual bits for bit-level cubes/clauses
**        However, the transition system itself still uses bitvectors
**/

#include "engines/ic3bits.h"

using namespace smt;
using namespace std;

namespace pono {

IC3Bits::IC3Bits(const SafetyProperty & p,
                 const TransitionSystem & ts,
                 const SmtSolver & s,
                 PonoOptions opt)
    : super(p, ts, s, opt)
{
}

void IC3Bits::initialize()
{
  if (super::initialized_) {
    return;
  }

  super::initialize();

  Term bv1 = solver_->make_term(1, solver_->make_sort(BV, 1));

  assert(!state_bits_.size());
  for (const auto & sv : ts_.statevars()) {
    const Sort & sort = sv->get_sort();
    if (sort == boolsort_) {
      state_bits_.push_back(sv);
    } else {
      assert(sort->get_sort_kind() == BV);
      for (size_t i = 0; i < sort->get_width(); ++i) {
        state_bits_.push_back(solver_->make_term(
            Equal, solver_->make_term(Op(Extract, i, i), sv), bv1));
      }
    }
  }
}

IC3Formula IC3Bits::get_model_ic3formula() const
{
  // expecting all solving in IC3 to be done at context level > 0
  // so if we're getting a model we should not be at context 0
  assert(solver_context_);

  TermVec children;
  children.reserve(state_bits_.size());
  for (const auto & b : state_bits_) {
    if (solver_->get_value(b) == solver_true_) {
      children.push_back(b);
    } else {
      children.push_back(solver_->make_term(Not, b));
    }
  }

  return ic3formula_conjunction(children);
}

bool IC3Bits::ic3formula_check_valid(const IC3Formula & u) const
{
  // check that children are booleans
  // with only a single variable
  UnorderedTermSet free_vars;
  Op op;
  for (const auto & c : u.children) {
    free_vars.clear();
    get_free_symbolic_consts(c, free_vars);
    if (c->get_sort() != boolsort_ || free_vars.size() > 1) {
      return false;
    }
  }

  // got through all checks without failing
  return true;
}

void IC3Bits::check_ts() const
{
  for (const auto & sv : ts_.statevars()) {
    const Sort & sort = sv->get_sort();
    if (sort != boolsort_ && sort->get_sort_kind() != BV) {
      throw PonoException("Unsupported variable sort in IC3Bits: "
                          + sv->to_string() + ":" + sort->to_string());
    }
  }

  for (const auto & iv : ts_.inputvars()) {
    const Sort & sort = iv->get_sort();
    if (sort != boolsort_ && sort->get_sort_kind() != BV) {
      throw PonoException("Unsupported variable sort in IC3Bits: "
                          + iv->to_string() + ":" + sort->to_string());
    }
  }
}

}  // namespace pono
