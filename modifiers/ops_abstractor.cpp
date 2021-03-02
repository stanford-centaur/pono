/*********************                                                  */
/*! \file ops_abstractor.cpp
** \verbatim
** Top contributors (to current version):
**   Ahmed Irfan, Makai Mann
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

OpsAbstractor::OpsAbstractor(const TransitionSystem & conc_ts,
                             TransitionSystem & abs_ts)
    : super(conc_ts, abs_ts),
      solver_(abs_ts_.solver()),
      abs_walker_(*this, &abstraction_cache_),
      conc_walker_(*this, &concretization_cache_),
      min_bw_(0)
{
}

Term OpsAbstractor::abstract(Term & t)
{
  return abs_walker_.visit(t);
}

Term OpsAbstractor::concrete(Term & t) { return conc_walker_.visit(t); }

void OpsAbstractor::set_ops_to_abstract(const UnorderedOpSet & ops_to_abstract)
{
  if (ops_to_abstract_.size() > 0) {
    throw PonoException(
        "OpsAbstractor::set_ops_to_abstract "
        "Cannot reset the set of operators to abstract");
  }

  ops_to_abstract_.insert(ops_to_abstract.begin(), ops_to_abstract.end());
}

void OpsAbstractor::do_abstraction()
{
  if (ops_to_abstract_.size() == 0) {
    throw PonoException(
        "OpsAbstractor::do_abstraction "
        "No operators to abstract."
        "Set them by calling set_ops_to_abstract");
  }

  if (!abs_ts_.is_functional()) {
    abs_ts_ = conc_ts_;
    RelationalTransitionSystem & abs_rts =
      static_cast<RelationalTransitionSystem &>(abs_ts_);

    Term init = conc_ts_.init();
    Term trans = conc_ts_.trans();
    Term abs_init = abstract(init);
    Term abs_trans = abstract(trans);

    abs_rts.set_init(abs_init);
    abs_rts.set_trans(abs_trans);
  } else {
    FunctionalTransitionSystem & abs_fts =
      static_cast<FunctionalTransitionSystem &>(abs_ts_);

    // add inputs and statevars
    for (const auto & v : conc_ts_.inputvars()) {
      abs_fts.add_inputvar(v);
    }
    for (const auto &v : conc_ts_.statevars()) {
      abs_fts.add_statevar(v, conc_ts_.next(v));
    }

    // state updates
    for (const auto & e : conc_ts_.state_updates()) {
      Term val = e.second;
      abs_fts.assign_next(e.first, abstract(val));
    }

    Term init = conc_ts_.init();
    abs_fts.set_init(abstract(init));

    for (const auto & e : conc_ts_.named_terms()) {
      abs_fts.name_term(e.first, e.second);
    }
    for (const auto &c : conc_ts_.constraints()) {
      Term cc = c.first;
      abs_fts.add_constraint(abstract(cc), c.second);
    }
  }
}

OpsAbstractor::AbstractionWalker::AbstractionWalker(
    OpsAbstractor & oa, UnorderedTermMap * ext_cache)
    : IdentityWalker(oa.solver_, false, ext_cache), oa_(oa)
{
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
  // check if we do not need to abstract the operator
  if (op.is_null() ||
      oa_.ops_to_abstract_.find(op) == oa_.ops_to_abstract_.end() ||
      (sk == BV && sort->get_width() <= oa_.min_bw_)) {
    res = op.is_null() ? term : solver_->make_term(op, cached_children);
  } else {
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
        string op_str =
            "abs_" + op.to_string() + to_string(sk) + "_"
            + to_string(cached_children[0]->get_sort()->get_sort_kind());

        SortVec sv({ cached_children[0]->get_sort() });
        if (cached_children.size() == 2) {
          op_str +=
              "_" + to_string(cached_children[1]->get_sort()->get_sort_kind());
          sv.push_back(cached_children[1]->get_sort());
        }
        sv.push_back(sort);

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

        TermVec args = { abs_op };
        args.insert(args.end(), cached_children.begin(), cached_children.end());
        res = solver_->make_term(Apply, args);

        assert(oa_.abs_terms_.find(res) == oa_.abs_terms_.end());
        oa_.abs_terms_[res] = solver_->make_term(op, cached_children);
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
        string op_str =
            "abs_" + op.to_string() + "_" + to_string(sort->get_width()) + "_"
            + to_string(cached_children[0]->get_sort()->get_width());

        SortVec sv({ cached_children[0]->get_sort() });
        if (cached_children.size() == 2) {
          op_str +=
              "_" + to_string(cached_children[1]->get_sort()->get_width());
          sv.push_back(cached_children[1]->get_sort());
        }
        sv.push_back(sort);

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

        TermVec args = { abs_op };
        args.insert(args.end(), cached_children.begin(), cached_children.end());
        res = solver_->make_term(Apply, args);

        assert(oa_.abs_terms_.find(res) == oa_.abs_terms_.end());
        oa_.abs_terms_[res] = solver_->make_term(op, cached_children);
        break;
      }
        // other operators
      default:
        // should not reach
        assert(false);
        break;
    };
  }

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
