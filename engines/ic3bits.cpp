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

#include <cassert>
#include <cstddef>

#include "engines/ic3base.h"
#include "smt-switch/smt.h"
#include "smt-switch/utils.h"
#include "utils/exceptions.h"
#include "utils/partial_model.h"

using namespace smt;

namespace pono {

IC3Bits::IC3Bits(const SafetyProperty & p,
                 const TransitionSystem & ts,
                 const SmtSolver & solver,
                 PonoOptions opt,
                 Engine engine)
    : super(p, ts, solver, opt, engine),
      partial_model_getter_(solver_),
      has_assumptions_(false)
{
}

void IC3Bits::initialize()
{
  if (super::initialized_) {
    return;
  }

  super::initialize();

  build_ts_related_info();

  // Constraints live over current/no-next terms. Include each constraint and
  // its next-state update rewritten back to current variables when computing
  // the predecessor's relevant bit slices.
  assert(!nxt_state_updates_.empty());
  for (const auto & c_initnext : ts_.constraints()) {
    has_assumptions_ = true;
    assert(ts_.no_next(c_initnext.first));
    constraints_curr_var_.emplace(
        c_initnext.first, std::vector<std::pair<int, int>>({ { 0, 0 } }));
    constraints_curr_var_.emplace(
        next_curr_replace(ts_.next(c_initnext.first)),
        std::vector<std::pair<int, int>>({ { 0, 0 } }));
  }

  Term bv1 = solver_->make_term(1, solver_->make_sort(BV, 1));

  assert(!state_bits_.size());
  for (const auto & sv : ts_.statevars()) {
    const Sort & sort = sv->get_sort();
    if (sort == boolsort_) {
      state_bits_.push_back(sv);
    } else {
      assert(sort->get_sort_kind() == BV);
      for (std::size_t i = 0; i < sort->get_width(); ++i) {
        state_bits_.push_back(solver_->make_term(
            Equal, solver_->make_term(Op(Extract, i, i), sv), bv1));
      }
    }
  }
}

void IC3Bits::build_ts_related_info()
{
  // nxt_state_updates_ rewrites primed state variables to their next-state
  // functions. State variables without updates are treated like inputs and are
  // omitted from partial predecessors unless assumptions force keeping them.
  const auto & state_updates = ts_.state_updates();
  for (const auto & sv : ts_.statevars()) {
    auto pos = state_updates.find(sv);
    if (pos == state_updates.end()) {
      no_next_vars_.insert(sv);
    } else {
      nxt_state_updates_.emplace(ts_.next(sv), pos->second);
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
  for (const auto & vars : { ts_.statevars(), ts_.inputvars() }) {
    for (const auto & var : vars) {
      const Sort & sort = var->get_sort();
      if (sort != boolsort_ && sort->get_sort_kind() != BV) {
        throw PonoException("IC3Bits only supports bit-vectors, got "
                            + to_string(sort->get_sort_kind()));
      }
    }
  }
}

bool IC3Bits::keep_var_in_partial_model(const Term & v) const
{
  // With constraints, input/no-next variables can affect enabled transitions,
  // so keep all current variables. Otherwise restrict cubes to updated state.
  if (has_assumptions_) {
    return ts_.is_curr_var(v);
  }

  return ts_.is_curr_var(v) && no_next_vars_.find(v) == no_next_vars_.end();
}

IC3Formula IC3Bits::ExtractPartialModel(const Term & p)
{
  assert(ts_.no_next(p));

  std::unordered_map<Term, std::vector<std::pair<int, int>>> varlist;
  // IC3 gives a current-state bad predicate. To compute a predecessor cube,
  // move it through the next-state functions and then analyze the result over
  // current variables.
  Term bad_state_no_nxt = next_curr_replace(ts_.next(p));

  if (has_assumptions_) {
    auto ast_slices = constraints_curr_var_;
    ast_slices.emplace(bad_state_no_nxt,
                       std::vector<std::pair<int, int>>({ { 0, 0 } }));
    partial_model_getter_.GetVarListForAsts_in_bitlevel(ast_slices, varlist);
  } else {
    std::unordered_map<Term, std::vector<std::pair<int, int>>> in_ast;
    in_ast.emplace(bad_state_no_nxt,
                   std::vector<std::pair<int, int>>({ { 0, 0 } }));
    partial_model_getter_.GetVarListForAsts_in_bitlevel(in_ast, varlist);
  }

  TermVec conjvec_partial;
  for (const auto & v_slice_pair : varlist) {
    const auto & v = v_slice_pair.first;
    if (!keep_var_in_partial_model(v)) continue;

    Term val = solver_->get_value(v);
    if (v->get_sort() == boolsort_) {
      // Boolean variables are already single-bit terms; Extract is only valid
      // for bit-vectors.
      conjvec_partial.push_back(solver_->make_term(Op(Equal), v, val));
      continue;
    }

    for (const auto & h_l_pair : v_slice_pair.second) {
      for (int idx = h_l_pair.first; idx >= h_l_pair.second; --idx) {
        auto eq = solver_->make_term(
            Op(Equal),
            solver_->make_term(Op(Extract, idx, idx), v),
            solver_->make_term(Op(Extract, idx, idx), val));
        conjvec_partial.push_back(eq);
      }
    }
  }

  if (conjvec_partial.empty()) conjvec_partial.push_back(solver_true_);
  return ic3formula_conjunction(conjvec_partial);
}

void IC3Bits::predecessor_generalization(size_t i,
                                         const Term & cterm,
                                         IC3Formula & pred)
{
  if (!options_.ic3_pregen_) return;

  if (options_.ic3bits_partial_model_pregen_) {
    pred = ExtractPartialModel(cterm);
    return;
  }

  IC3::predecessor_generalization(i, cterm, pred);
}

}  // namespace pono
