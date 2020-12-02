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

#include <random>

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

IC3Unit ClauseHandler::create_negated(const smt::TermVec & c) const
{
  assert(c.size());
  Term term = c.at(0);
  for (size_t i = 1; i < c.size(); ++i) {
    term = solver_->make_term(And, term, c[i]);
  }
  IC3Unit res(term, c, true);
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

/** IC3 Implementation */

IC3::IC3(Property & p, smt::SolverEnum se)
    : super(p, se, unique_ptr<ClauseHandler>(new ClauseHandler(solver_)))
{
}

IC3::IC3(Property & p, const smt::SmtSolver & s)
    : super(p, s, unique_ptr<ClauseHandler>(new ClauseHandler(solver_)))
{
}

IC3::IC3(const PonoOptions & opt, Property & p, smt::SolverEnum se)
    : super(opt, p, se, unique_ptr<ClauseHandler>(new ClauseHandler(solver_)))
{
}

IC3::IC3(const PonoOptions & opt, Property & p, const smt::SmtSolver & s)
    : super(opt, p, s, unique_ptr<ClauseHandler>(new ClauseHandler(solver_)))
{
}

std::vector<IC3Unit> IC3::inductive_generalization(size_t i, const IC3Unit & c)
{
  assert(c.negated);  // expecting a cube

  if (options_.ic3_indgen_mode_ != 0) {
    throw PonoException("Boolean IC3 only supports indgen mode 0 but got "
                        + std::to_string(options_.ic3_indgen_mode_));
  }

  UnorderedTermSet keep, core_set;
  TermVec bool_assump, tmp, new_tmp, removed, lits;
  lits = c.children;

  if (options_.random_seed_ > 0) {
    shuffle(
        lits.begin(), lits.end(), default_random_engine(options_.random_seed_));
  }

  int iter = 0;
  bool progress = true;
  while (iter <= options_.ic3_gen_max_iter_ && lits.size() > 1 && progress) {
    iter = options_.ic3_gen_max_iter_ > 0 ? iter + 1 : iter;
    size_t prev_size = lits.size();
    for (auto a : lits) {
      // check if we can drop a
      if (keep.find(a) != keep.end()) {
        continue;
      }
      tmp.clear();
      for (auto aa : lits) {
        if (a != aa) {
          tmp.push_back(aa);
        }
      }

      Term tmp_and_term = make_and(tmp);
      if (!intersects_initial(tmp_and_term)) {
        assert(solver_context_ == 0);
        push_solver_context();
        assert_frame_labels(i - 1);
        assert_trans_label();
        solver_->assert_formula(solver_->make_term(Not, tmp_and_term));

        Term l;
        bool_assump.clear();
        for (auto t : tmp) {
          l = label(t);
          solver_->assert_formula(solver_->make_term(Implies, l, ts_.next(t)));
          bool_assump.push_back(l);
        }

        Result r = solver_->check_sat_assuming(bool_assump);
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

  IC3Unit blocking_clause = handler_->create(lits);
  assert(!blocking_clause.negated);  // expecting a clause
  return { blocking_clause };
}

IC3Unit IC3::generalize_predecessor(size_t i, const IC3Unit & c)
{
  assert(i > 0);
  if (i == 1) {
    // don't need to generalize if i == 1
    // the predecessor is an initial state
    assert(intersects_initial(c.term));
    return c;
  }

  // should never intersect with a frame before F[i-1]
  // otherwise, this predecessor should have been found
  // in a previous step (before a new frame was pushed)
  assert(i > 1);  // so F[i-2] makes sense
  assert(!intersects(c.term, get_frame(i - 2)));

  const UnorderedTermSet & statevars = ts_.statevars();
  TermVec next_lits;
  next_lits.reserve(statevars.size());
  for (auto v : statevars) {
    Term nv = ts_.next(v);
    assert(ts_.is_next_var(nv));
    Term next_val = solver_->get_value(nv);
    next_lits.push_back(solver_->make_term(Equal, nv, next_val));
  }

  // collect input assignments
  const UnorderedTermSet & inputvars = ts_.inputvars();
  TermVec input_lits;
  input_lits.reserve(inputvars.size());
  for (auto v : inputvars) {
    Term val = solver_->get_value(v);
    input_lits.push_back(solver_->make_term(Equal, v, val));
  }

  Term formula = make_and(input_lits);
  if (ts_.is_deterministic()) {
    formula = solver_->make_term(And, formula, trans_label_);
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
    Term pre_formula = get_frame(i - 1);
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
      formula, c.children, red_cube_lits, &rem_cube_lits);

  // should need some assumptions
  // formula should not be unsat on its own
  assert(red_cube_lits.size() > 0);

  IC3Unit res = handler_->create_negated(red_cube_lits);
  // expecting a Cube here
  assert(res.negated);

  return res;
}

void IC3::check_ts() const
{
  Sort boolsort = solver_->make_sort(BOOL);
  for (auto sv : ts_.statevars()) {
    if (sv->get_sort() != boolsort) {
      throw PonoException("Got non-boolean state variable in bit-level IC3: "
                          + sv->to_string());
    }
  }

  for (auto iv : ts_.inputvars()) {
    if (iv->get_sort() != boolsort) {
      throw PonoException("Got non-boolean input variable in bit-level IC3: "
                          + iv->to_string());
    }
  }
}

IC3Unit IC3::get_unit() const
{
  // expecting all solving in IC3 to be done at context level > 0
  // so if we're getting a model we should not be at context 0
  assert(solver_context_);

  const UnorderedTermSet & statevars = ts_.statevars();
  TermVec children;
  children.reserve(statevars.size());
  for (auto sv : ts_.statevars()) {
    if (solver_->get_value(sv) == solver_true_) {
      children.push_back(sv);
    } else {
      children.push_back(solver_->make_term(Not, sv));
    }
  }

  IC3Unit res = handler_->create_negated(children);
  assert(res.negated);  // expecting a Cube here
  return res;
}

}  // namespace pono
