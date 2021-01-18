/*********************                                                  */
/*! \file ops_abstractor.cpp
** \verbatim
** Top contributors (to current version):
**   Ahmed Irfan
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Abstract arrays using uninterpreted functions.
**
**
**/

#include "modifiers/ops_abstractor.h"

using namespace smt;
using namespace std;

namespace pono {

OpsAbstractor::OpsAbstractor(const TransitionSystem & conc_ts, TransitionSystem & abs_ts)
  : super(conc_ts, abs_ts),
    solver_(abs_ts_.solver()),
    abs_walker_(*this, &abstraction_cache_),
    conc_walker_(*this, &concretization_cache_)
{
}

Term OpsAbstractor::abstract(Term & t)
{
  return abs_walker_.visit(t);
}

Term  OpsAbstractor::concrete(Term & t)
{
  return conc_walker_.visit(t);
}

OpsAbstractor::AbstractionWalker::AbstractionWalker(OpsAbstractor &oa,
                                                    UnorderedTermMap *ext_cache)
  : IdentityWalker(oa.solver_, false, ext_cache),
    oa_(oa)
{
}

void OpsAbstractor::do_abstraction()
{
  Term init = conc_ts_.init();
  Term trans = conc_ts_.trans();
  Term abs_init = abstract(init);
  Term abs_trans = abstract(trans);

  // for now, supporting a relational system
  // but generic abstractor does not require a relational system
  // do a cast
  assert(!abs_ts_.is_functional());
  RelationalTransitionSystem & abs_rts =
    static_cast<RelationalTransitionSystem &>(abs_ts_);
  // the calls to abstract have already added the
  // (possibly abstracted) variables
  // now we just need to set the initial states and trans
  abs_rts.set_init(abs_init);
  abs_rts.set_trans(abs_trans);
}

WalkerStepResult OpsAbstractor::AbstractionWalker::visit_term(Term & term)
{
  if (preorder_ || in_cache(term)) {
    return Walker_Continue;
  }

  Sort sort = term->get_sort();
  SortKind sk = sort->get_sort_kind();
  Op op = term->get_op();

  TermVec cached_children;
  Term cc;
  for (auto c : term) {
    bool ok = query_cache(c, cc);
    assert(ok);  // in post-order so should always have a cache hit
    cached_children.push_back(cc);
  }

  Term res;
  switch (op.prim_op) {
    // arithmetic operators
  case Plus:
  case Minus:
  case Mult:
  case Div:
  case Mod:
  case Pow:
  case IntDiv:
  case Abs: {
    assert(cached_children.size() <= 2);
    string op_str = "abs_" + op.to_string() + to_string(sk) + "_"
      + to_string(cached_children[0]->get_sort()->get_sort_kind());
    SortVec sv({sort, cached_children[0]->get_sort()});
    if (cached_children.size() == 2) {
      op_str += "_" + to_string(cached_children[1]->get_sort()->get_sort_kind());
      sv.push_back(cached_children[1]->get_sort());
    }

    Term abs_op;
    auto it = oa_.abs_op_symbols_.find(op_str);
    if (it == oa_.abs_op_symbols_.end()) {
      Sort func_sort = solver_->make_sort(FUNCTION, sv);
      abs_op = solver_->make_symbol(op_str, func_sort);
      oa_.abs_op_symbols_[op_str] = abs_op;
      oa_.abs_symbols_to_op_[abs_op] = op;
    } else {
      abs_op = it->second;
    }

    TermVec args = {abs_op};
    args.insert(args.end(), cached_children.begin(), cached_children.end());
    res = solver_->make_term(Apply, args);
    break;
  }
    // BV operators
  case Concat:
  case Extract:
  case BVAnd:
  case BVOr:
  case BVXor:
  case BVNand:
  case BVNor:
  case BVXnor:
  case BVAdd:
  case BVSub:
  case BVMul:
  case BVUdiv:
  case BVSdiv:
  case BVUrem:
  case BVSrem:
  case BVSmod:
  case BVShl:
  case BVAshr:
  case BVLshr:
  case BVNot:
  case BVNeg: {
    assert(cached_children.size() <= 2);
    string op_str = "abs_" + op.to_string() + "_" + to_string(sort->get_width()) + "_"
      + to_string(cached_children[0]->get_sort()->get_width());
    SortVec sv({sort, cached_children[0]->get_sort()});
    if (cached_children.size() == 2) {
      op_str += "_" + to_string(cached_children[1]->get_sort()->get_width());
      sv.push_back(cached_children[1]->get_sort());
    }

    Term abs_op;
    auto it = oa_.abs_op_symbols_.find(op_str);
    if (it == oa_.abs_op_symbols_.end()) {
      Sort func_sort = solver_->make_sort(FUNCTION, sv);
      abs_op = solver_->make_symbol(op_str, func_sort);
      oa_.abs_op_symbols_[op_str] = abs_op;
      oa_.abs_symbols_to_op_[abs_op] = op;
    } else {
      abs_op = it->second;
    }

    TermVec args = {abs_op};
    args.insert(args.end(), cached_children.begin(), cached_children.end());
    res = solver_->make_term(Apply, args);
    break;
  }
    // other operators
  default:
    if (op.is_null()) {
      res = term;
    } else {
      res = solver_->make_term(op, cached_children);
    }
    break;
  };

  assert(res);
  oa_.update_term_cache(term, res);

  return Walker_Continue;
}

OpsAbstractor::ConcretizationWalker::ConcretizationWalker(OpsAbstractor &oa,
                                                       UnorderedTermMap *ext_cache)
  : IdentityWalker(oa.solver_, false, ext_cache),
    oa_(oa)
{
}

WalkerStepResult OpsAbstractor::ConcretizationWalker::visit_term(Term & term)
{
  if (preorder_ || in_cache(term)) {
    return Walker_Continue;
  }

  Op op = term->get_op();

  // if not an apply, then we don't need to do anything except
  // rebuild to keep changes
  if (op != Apply) {
    if (op.is_null()) {
      oa_.update_term_cache(term, term);
    } else {
      TermVec cached_children;
      Term cc;
      for (auto c : term) {
        bool ok = query_cache(c, cc);
        assert(ok);
        cached_children.push_back(cc);
      }

      Term rebuilt = solver_->make_term(op, cached_children);
      // Note: reversed order because update_term_cache takes the
      // concrete term first
      oa_.update_term_cache(rebuilt, term);
    }
    return Walker_Continue;
  }

  assert(op == Apply);
  auto it = term->begin();
  Term uf = *it;
  TermVec cached_args;
  ++it;
  while (it != term->end()) {
    Term ca;
    bool ok = query_cache(*it, ca);
    assert(ok);
    assert(ca);
    cached_args.push_back(ca);
    ++it;
  }

  Term res;
  if (oa_.abs_symbols_to_op_.find(uf) != oa_.abs_symbols_to_op_.end()) {
    res = solver_->make_term(oa_.abs_symbols_to_op_.at(uf), cached_args);
  } else {
    cached_args.insert(cached_args.begin(), uf);
    res = solver_->make_term(op, cached_args);
  }

  assert(res);
  oa_.update_term_cache(res, term);

  return Walker_Continue;
}

} // namespace pono
