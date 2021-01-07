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

#include <algorithm>
#include <random>

#include "assert.h"
#include "utils/term_analysis.h"

using namespace smt;
using namespace std;

namespace pono {

/** IC3 Implementation */

IC3::IC3(const Property & p, const TransitionSystem & ts,
         const smt::SmtSolver & s, PonoOptions opt)
  : super(p, ts, s, opt)
{
  engine_ = Engine::IC3_BOOL;
  solver_->set_opt("produce-unsat-cores", "true");
}

IC3Formula IC3::get_model_ic3formula(TermVec * out_inputs,
                                     TermVec * out_nexts) const
{
  // expecting all solving in IC3 to be done at context level > 0
  // so if we're getting a model we should not be at context 0
  assert(solver_context_);

  const UnorderedTermSet & statevars = ts_.statevars();
  TermVec children;
  children.reserve(statevars.size());
  for (const auto &sv : ts_.statevars()) {
    if (solver_->get_value(sv) == solver_true_) {
      children.push_back(sv);
    } else {
      children.push_back(solver_->make_term(Not, sv));
    }

    if (out_nexts) {
      Term nv = ts_.next(sv);
      if (solver_->get_value(nv) == solver_true_) {
        out_nexts->push_back(nv);
      } else {
        out_nexts->push_back(solver_->make_term(Not, nv));
      }
    }
  }

  if (out_inputs) {
    for (const auto &iv : ts_.inputvars()) {
      if (solver_->get_value(iv) == solver_true_) {
        out_inputs->push_back(iv);
      } else {
        out_inputs->push_back(solver_->make_term(Not, iv));
      }
    }
  }

  return ic3formula_conjunction(children);
}

bool IC3::ic3formula_check_valid(const IC3Formula & u) const
{
  const Sort &boolsort = solver_->make_sort(BOOL);
  // check that children are literals
  Op op;
  for (const auto &c : u.children) {
    if (!is_lit(c, boolsort)) {
      return false;
    }
  }

  // for now not checking the actual term (e.g. u.term)
  // it's somewhat hard if the underlying solver does rewriting

  // got through all checks without failing
  return true;
}

IC3Formula IC3::inductive_generalization(size_t i, const IC3Formula & c)
{
  assert(!c.disjunction);  // expecting a cube

  if (options_.mbic3_indgen_mode != 0) {
    throw PonoException("Boolean IC3 only supports indgen mode 0 but got "
                        + std::to_string(options_.mbic3_indgen_mode));
  }

  UnorderedTermSet keep, core_set;
  TermVec bool_assump, tmp, new_tmp, removed, lits;
  lits = c.children;

  if (options_.random_seed_ > 0) {
    shuffle(
        lits.begin(), lits.end(), default_random_engine(options_.random_seed_));
  }

  //TODO: use unsatcore-reducer
  int iter = 0;
  bool progress = true;
  while (iter <= options_.ic3_gen_max_iter_ && lits.size() > 1 && progress) {
    iter = options_.ic3_gen_max_iter_ > 0 ? iter + 1 : iter;
    size_t prev_size = lits.size();
    for (const auto &a : lits) {
      // check if we can drop a
      if (keep.find(a) != keep.end()) {
        continue;
      }
      tmp.clear();
      for (const auto &aa : lits) {
        if (a != aa) {
          tmp.push_back(aa);
        }
      }

      Term tmp_and_term = make_and(tmp);
      if (!check_intersects_initial(tmp_and_term)) {
        assert(solver_context_ == 0);
        push_solver_context();
        assert_frame_labels(i - 1);
        assert_trans_label();
        solver_->assert_formula(solver_->make_term(Not, tmp_and_term));

        Term l;
        bool_assump.clear();
        for (const auto &t : tmp) {
          l = label(t);
          solver_->assert_formula(solver_->make_term(Implies, l, ts_.next(t)));
          bool_assump.push_back(l);
        }

        Result r = check_sat_assuming(bool_assump);
        assert(!r.is_unknown());

        if (r.is_sat()) {
          // we cannot drop a
          pop_solver_context();
        } else {
          new_tmp.clear();
          removed.clear();
          core_set.clear();
          // filter using unsatcore
          solver_->get_unsat_core(core_set);
          for (size_t j = 0; j < bool_assump.size(); ++j) {
            if (core_set.find(bool_assump[j]) != core_set.end()) {
              new_tmp.push_back(tmp[j]);
            } else {
              removed.push_back(tmp[j]);
            }
          }

          pop_solver_context();

          // keep in mind that you cannot drop a literal if it causes c to
          // intersect with the initial states
          size_t size = new_tmp.size();
          fix_if_intersects_initial(new_tmp, removed);
          // remember the literals which cannot be dropped
          for (size_t i = size; i < new_tmp.size(); ++i) {
            keep.insert(new_tmp[i]);
          }

          lits = new_tmp;
          break;  // next iteration
        }
      }
    }

    progress = lits.size() < prev_size;
  }

  // TODO: would it be more intuitive to start with a clause
  //       and generalize the clause directly?
  const IC3Formula &blocking_clause = ic3formula_negate(ic3formula_conjunction(lits));
  assert(blocking_clause.disjunction);  // expecting a clause
  return blocking_clause;
}

IC3Formula IC3::generalize_predecessor(size_t i, const IC3Formula & c)
{
  // TODO: change this so we don't have to depend on the solver context to be
  // sat
  assert(i > 0);

  const UnorderedTermSet & statevars = ts_.statevars();
  TermVec input_lits, next_lits;
  const IC3Formula & icf = get_model_ic3formula(&input_lits, &next_lits);
  const TermVec & cube_lits = icf.children;

  if (i == 1) {
    // don't need to generalize if i == 1
    // the predecessor is an initial state
    return get_model_ic3formula();
  }

  Term formula = make_and(input_lits);
  if (ts_.is_deterministic()) {
    // NOTE: need to use full trans, not just trans_label_ here
    //       because we are passing it to the reducer_
    formula = solver_->make_term(And, formula, ts_.trans());
    formula = solver_->make_term(
        And, formula, solver_->make_term(Not, ts_.next(c.term)));
  } else {
    formula = solver_->make_term(And, formula, make_and(next_lits));

    // preimage formula
    // NOTE: because this will be negated, it is important to use the
    // whole frame and trans, not just the labels
    // because label is: trans_label_ -> trans
    // so if it is negated, that doesn't force trans to be false
    // the implication could be more efficient than iff so we want to leave it
    // that way
    Term pre_formula = get_frame_term(i - 1);
    pre_formula = solver_->make_term(And, pre_formula, ts_.trans());
    pre_formula =
        solver_->make_term(And, pre_formula, solver_->make_term(Not, c.term));
    pre_formula = solver_->make_term(And, pre_formula, ts_.next(c.term));
    formula =
        solver_->make_term(And, formula, solver_->make_term(Not, pre_formula));
  }
  // TODO: consider adding functional pre-image here
  //       not sure if it makes sense to have at the boolean level

  TermVec red_cube_lits, rem_cube_lits;
  reducer_.reduce_assump_unsatcore(
      formula, cube_lits, red_cube_lits, &rem_cube_lits);

  // should need some assumptions
  // formula should not be unsat on its own
  assert(red_cube_lits.size() > 0);

  const IC3Formula & res = ic3formula_conjunction(red_cube_lits);
  // expecting a Cube here
  assert(!res.disjunction);

  return res;
}

void IC3::check_ts() const
{
  const Sort &boolsort = solver_->make_sort(BOOL);
  for (const auto &sv : ts_.statevars()) {
    if (sv->get_sort() != boolsort) {
      throw PonoException("Got non-boolean state variable in bit-level IC3: "
                          + sv->to_string());
    }
  }

  for (const auto &iv : ts_.inputvars()) {
    if (iv->get_sort() != boolsort) {
      throw PonoException("Got non-boolean input variable in bit-level IC3: "
                          + iv->to_string());
    }
  }
}

}  // namespace pono
