/*********************                                                  */
/*! \file mbic3.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Simple implementation of IC3 using model values
**
**/

#include "engines/mbic3.h"

#include "smt-switch/utils.h"
#include "smt/available_solvers.h"
#include "utils/logger.h"

using namespace smt;
using namespace std;

namespace pono {

static const std::unordered_set<smt::PrimOp> expected_ops(
    { Equal, BVComp, Ge, Gt, Le, Lt, BVUle, BVUlt, BVUge, BVUgt });

static bool disjoint_set_rank(const Term & t1, const Term & t2)
{
  if (!t1->is_value() && !t2->is_value()) {
    return t1 < t2;
  }
  return !t1->is_value();
}

static void split_eq(SmtSolver & solver, const TermVec & in, TermVec & out)
{
  TermVec args;
  for (auto a : in) {
    if (a->get_op().prim_op == PrimOp::Equal) {
      args.clear();
      for (auto aa : a) {
        args.push_back(aa);
      }

      Sort s = args[0]->get_sort();
      SortKind sk = s->get_sort_kind();
      if (sk == SortKind::BOOL) {
        out.push_back(a);  // assert(false);
      } else if (sk == SortKind::BV) {
        out.push_back(solver->make_term(BVUle, args[0], args[1]));
        out.push_back(solver->make_term(BVUle, args[1], args[0]));
      } else if (sk == SortKind::INT || sk == SortKind::REAL) {
        out.push_back(solver->make_term(Le, args[0], args[1]));
        out.push_back(solver->make_term(Le, args[1], args[0]));
      } else {
        out.push_back(a);
      }
    } else {
      out.push_back(a);
    }
  }
}

ModelBasedIC3::ModelBasedIC3(const SafetyProperty & p,
                             const TransitionSystem & ts,
                             const SmtSolver & slv,
                             PonoOptions opt)
    : super(p, ts, slv, opt)
{
  engine_ = Engine::MBIC3;
  interpolator_enum_ = opt.smt_interpolator_;
}

IC3Formula ModelBasedIC3::get_model_ic3formula() const
{
  DisjointSet ds(disjoint_set_rank);
  TermVec cube_lits;
  const UnorderedTermSet & statevars = ts_.statevars();

  for (const auto & v : statevars) {
    Term val = solver_->get_value(v);
    cube_lits.push_back(solver_->make_term(Equal, v, val));
    ds.add(v, val);
    assert(ts_.is_curr_var(v));
  }

  // add equalities from disjoint set
  for (const auto & v : statevars) {
    Term t = ds.find(v);
    if (t != v) {
      cube_lits.push_back(solver_->make_term(Equal, t, v));
    }
  }

  return ic3formula_conjunction(cube_lits);
}

bool ModelBasedIC3::ic3formula_check_valid(const IC3Formula & u) const
{
  const Sort & boolsort = solver_->make_sort(BOOL);
  // check that children are literals
  Op op;
  for (const auto & c : u.children) {
    if (c->get_sort() != boolsort) {
      return false;
    }
  }

  // not checking that they are equalities etc...
  // hard with rewriting
  // but more specifically, it's not always true
  // in this IC3 version, we add the whole bad as a proof goal
  // (or at least the conjuncts in a conjunctive partition of bad)
  // regardless of the shape

  // got through all checks without failing
  return true;
}

IC3Formula ModelBasedIC3::inductive_generalization(size_t i,
                                                   const IC3Formula & c)
{
  // only need interpolator if in mode 2
  assert(options_.mbic3_indgen_mode == 2 || interpolator_ == nullptr);

  assert(!c.disjunction);
  IC3Formula gen_res;

  if (options_.ic3_indgen_) {
    if (options_.mbic3_indgen_mode == 0) {
      return super::inductive_generalization(i, c);
    } else if (options_.mbic3_indgen_mode == 1) {
      TermVec tmp, lits, red_lits;
      for (const auto & a : c.children) {
        tmp.push_back(ts_.next(a));
      }
      split_eq(solver_, tmp, lits);

      // ( (frame /\ trans /\ not(c)) \/ init') /\ c' is unsat
      Term formula = make_and({ get_frame_term(i - 1),
                                ts_.trans(),
                                solver_->make_term(Not, c.term) });
      formula = solver_->make_term(Or, formula, ts_.next(ts_.init()));
      reducer_.reduce_assump_unsatcore(formula,
                                       lits,
                                       red_lits,
                                       NULL,
                                       options_.ic3_gen_max_iter_,
                                       options_.random_seed_);
      TermVec curr_lits;
      curr_lits.reserve(red_lits.size());
      for (const auto & l : red_lits) {
        curr_lits.push_back(ts_.curr(l));
      }
      gen_res = ic3formula_negate(ic3formula_conjunction(curr_lits));
    } else if (options_.mbic3_indgen_mode == 2) {
      // TODO: consider creating a separate derived class for this mode
      interpolator_->reset_assertions();

      TermVec conjuncts;
      split_eq(solver_, c.children, conjuncts);

      // ( (frame /\ trans /\ not(c)) \/ init') /\ c' is unsat
      Term formula = make_and({ get_frame_term(i - 1),
                                ts_.trans(),
                                solver_->make_term(Not, make_and(conjuncts)) });
      formula = solver_->make_term(Or, formula, ts_.next(ts_.init()));

      Term int_A = to_interpolator_->transfer_term(formula, BOOL);
      // still use c in B
      // only split equalities in A to encourage more general unsat proofs /
      // interpolants
      Term int_B = to_interpolator_->transfer_term(ts_.next(c.term), BOOL);

      Term interp;
      Result r = interpolator_->get_interpolant(int_A, int_B, interp);
      assert(r.is_unsat());
      logger.log(3, "Got interpolant: {}", interp);

      // Note: structure of interpolant is not guaranteed to be a clause
      //       but we'll at least try to split it up into disjuncts
      // Note 2: using disjunctive_partition before translating back to
      //         main solver because boolector removes ORs when putting into
      //         AIG form
      //         but since boolector doesn't support interpolation we know
      //         the interpolator is not boolector

      TermVec interp_disjuncts;
      disjunctive_partition(interp, interp_disjuncts);

      TermVec solver_interp_disjuncts;
      solver_interp_disjuncts.reserve(interp_disjuncts.size());

      Sort interp_boolsort = interpolator_->make_sort(BOOL);
      for (const auto & ii : interp_disjuncts) {
        assert(ii->get_sort() == interp_boolsort);
        solver_interp_disjuncts.push_back(
            ts_.curr(to_solver_->transfer_term(ii, BOOL)));
      }

      assert(solver_interp_disjuncts.size() == interp_disjuncts.size());
      gen_res = ic3formula_disjunction(solver_interp_disjuncts);
    } else {
      assert(false);
    }
  }
  assert(gen_res.term);
  return gen_res;
}

void ModelBasedIC3::predecessor_generalization(size_t i,
                                               const Term & c,
                                               IC3Formula & pred)
{
  // NOTE: for now this implementation doesn't use pred
  //       except to assign to it at the end
  //       need the model in a particular format
  assert(solver_context_ == 1);  // shouldn't use solver, solving all in
                                 // reducer_
  DisjointSet ds(disjoint_set_rank);
  UnorderedTermMap model;
  const UnorderedTermSet & statevars = ts_.statevars();

  TermVec cube_lits;
  cube_lits.reserve(statevars.size());
  TermVec next_lits;
  next_lits.reserve(statevars.size());

  for (const auto & v : statevars) {
    Term val = solver_->get_value(v);
    cube_lits.push_back(solver_->make_term(Equal, v, val));
    ds.add(v, val);
    assert(ts_.is_curr_var(v));
    assert(model.find(v) == model.end());
    model[v] = val;

    Term nv = ts_.next(v);
    assert(ts_.is_next_var(nv));
    Term next_val = solver_->get_value(nv);
    next_lits.push_back(solver_->make_term(Equal, nv, next_val));
    assert(model.find(nv) == model.end());
    model[nv] = next_val;
  }

  // collect input assignments
  const UnorderedTermSet & inputvars = ts_.inputvars();
  TermVec input_lits;
  input_lits.reserve(inputvars.size());
  for (const auto & v : inputvars) {
    Term val = solver_->get_value(v);
    input_lits.push_back(solver_->make_term(Equal, v, val));
    assert(model.find(v) == model.end());
    model[v] = val;
  }

  assert(i > 0);
  if (i == 1) {
    // don't need to generalize if i == 1
    // the predecessor is an initial state
    return;
  }

  if (options_.ic3_pregen_ && !options_.ic3_functional_preimage_) {
    // add congruent equalities to cube_lits
    for (const auto & v : statevars) {
      Term t = ds.find(v);
      if (t != v) {
        cube_lits.push_back(solver_->make_term(Equal, t, v));
      }
    }

    Term formula = make_and(input_lits);

    if (ts_.is_deterministic()) {
      // NOTE: reducer doesn't have semantics for trans_label_
      // better to just use whole trans for now
      // Optionally in the future we could add those semantics
      // to the reducer_'s solver
      formula = solver_->make_term(And, formula, ts_.trans());
      formula = solver_->make_term(
          And, formula, solver_->make_term(Not, ts_.next(c)));
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
          solver_->make_term(And, pre_formula, solver_->make_term(Not, c));
      pre_formula = solver_->make_term(And, pre_formula, ts_.next(c));

      formula = solver_->make_term(
          And, formula, solver_->make_term(Not, pre_formula));
    }

    TermVec splits, red_cube_lits, rem_cube_lits;
    split_eq(solver_, cube_lits, splits);
    reducer_.reduce_assump_unsatcore(formula,
                                     splits,
                                     red_cube_lits,
                                     &rem_cube_lits,
                                     options_.ic3_gen_max_iter_,
                                     options_.random_seed_);
    // should need some assumptions
    // formula should not be unsat on its own
    assert(red_cube_lits.size() > 0);

    // update pred to the generalization
    pred = ic3formula_conjunction(red_cube_lits);

  } else if (options_.ic3_pregen_ && options_.ic3_functional_preimage_) {
    assert(ts_.is_deterministic());

    UnorderedTermMap m;
    for (const auto & v : inputvars) {
      m[v] = model.at(v);
    }
    for (const auto & v : statevars) {
      Term nv = ts_.next(v);
      m[nv] = model.at(nv);
    }

    Term fun_preimage = solver_->substitute(trans_label_, m);
    TermVec conjuncts;
    conjunctive_partition(fun_preimage, conjuncts, true);
    pred = ic3formula_conjunction(conjuncts);
  }
}

void ModelBasedIC3::check_ts() const
{
  // check if there are arrays or uninterpreted sorts and fail if so
  for (const auto & vec : { ts_.statevars(), ts_.inputvars() }) {
    for (const auto & st : vec) {
      SortKind sk = st->get_sort()->get_sort_kind();
      if (sk == ARRAY) {
        throw PonoException("ModelBasedIC3 does not support arrays yet");
      } else if (sk == UNINTERPRETED) {
        throw PonoException(
            "ModelBasedIC3 does not support uninterpreted sorts yet.");
      }
    }
  }
}

void ModelBasedIC3::initialize()
{
  if (initialized_) {
    return;
  }

  super::initialize();

  // only need interpolator infrastructure for mode 2 (interpolation)
  if (options_.mbic3_indgen_mode == 2) {
    interpolator_ =
        create_interpolating_solver_for(interpolator_enum_, engine_);
    to_interpolator_ = std::make_unique<TermTranslator>(interpolator_);
    to_solver_ = std::make_unique<TermTranslator>(solver_);

    UnorderedTermMap & cache = to_solver_->get_cache();
    Term ns;
    for (const auto & s : ts_.statevars()) {
      // common variables are next states, unless used for refinement in IC3IA
      // then will refer to current state variables after untiming
      // need to cache both
      cache[to_interpolator_->transfer_term(s)] = s;
      ns = ts_.next(s);
      cache[to_interpolator_->transfer_term(ns)] = ns;
    }

    // need to add uninterpreted functions as well
    // first need to find them all
    // NOTE need to use get_free_symbols NOT get_free_symbolic_consts
    // because the latter ignores uninterpreted functions
    UnorderedTermSet free_symbols;
    get_free_symbols(ts_.init(), free_symbols);
    get_free_symbols(ts_.trans(), free_symbols);
    get_free_symbols(bad_, free_symbols);

    for (const auto & s : free_symbols) {
      assert(s->is_symbol());
      if (s->is_symbolic_const()) {
        // ignore constants
        continue;
      }
      cache[to_interpolator_->transfer_term(s)] = s;
    }
  }
}

}  // namespace pono
