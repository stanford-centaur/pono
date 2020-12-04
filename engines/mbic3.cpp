/*********************                                         */
/*! \file mbic3.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan, Florian Lonsing
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Simple implementation of IC3 using model values
**
**/

#include <algorithm>
#include <random>

#include "smt-switch/utils.h"

#include "engines/mbic3.h"
#include "smt/available_solvers.h"
#include "utils/logger.h"
#include "utils/term_analysis.h"
#include "utils/term_walkers.h"

using namespace smt;
using namespace std;

namespace pono {

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

// Conjunction implementations

/** Less than comparison of the hash of two terms
 *  for use in sorting
 *  @param t0 the first term
 *  @param t1 the second term
 *  @return true iff t0's hash is less than t1's hash
 */
bool term_hash_lt(const smt::Term & t0, const smt::Term & t1)
{
  return (t0->hash() < t1->hash());
}

Conjunction::Conjunction(const smt::SmtSolver & solver,
                         const smt::TermVec & conjuncts)
    : conjuncts_(conjuncts)
{
  // sort literals
  std::sort(conjuncts_.begin(), conjuncts_.end(), term_hash_lt);

  // shouldn't have an empty cube
  assert(conjuncts_.size());

  // create term
  term_ = conjuncts_[0];
  for (size_t i = 1; i < conjuncts_.size(); ++i) {
    term_ = solver->make_term(smt::And, term_, conjuncts_[i]);
  }
}

// main implementations

ModelBasedIC3::ModelBasedIC3(Property & p, SolverEnum se)
    : super(p, se),
      solver_true_(solver_->make_term(true)),
      solver_false_(solver_->make_term(false)),
      solver_context_(0)
{
  initialize();
}

ModelBasedIC3::ModelBasedIC3(Property & p, const SmtSolver & slv)
    : super(p, slv),
      solver_true_(solver_->make_term(true)),
      solver_false_(solver_->make_term(false)),
      solver_context_(0)
{
}

ModelBasedIC3::ModelBasedIC3(const PonoOptions & opt,
                             Property & p,
                             const SolverEnum se)
    : super(opt, p, se),
      solver_true_(solver_->make_term(true)),
      solver_false_(solver_->make_term(false)),
      solver_context_(0)
{
}

ModelBasedIC3::ModelBasedIC3(const PonoOptions & opt,
                             Property & p,
                             const SmtSolver & slv)
    : super(opt, p, slv),
      solver_true_(solver_->make_term(true)),
      solver_false_(solver_->make_term(false)),
      solver_context_(0)
{
}

ModelBasedIC3::~ModelBasedIC3() {}

void ModelBasedIC3::initialize()
{
  if (initialized_) {
    return;
  }

  super::initialize();

  // super sets other options
  solver_->set_opt("produce-unsat-cores", "true");

  frames_.clear();
  frame_labels_.clear();
  proof_goals_.clear();
  // first frame is always the initial states
  push_frame();
  constrain_frame(0, ts_.init());
  push_frame();

  assert(!interpolator_);
  assert(solver_);
  assert(!to_interpolator_);
  assert(!to_solver_);

  if (options_.ic3_indgen_mode_ == 2) {
    // TODO make an option to set the interpolator from outside
    initialize_interpolator();
  }
}

void ModelBasedIC3::initialize_interpolator()
{
  interpolator_ = create_interpolating_solver(MSAT_INTERPOLATOR);
  to_interpolator_ = std::make_unique<TermTranslator>(interpolator_);
  to_solver_ = std::make_unique<TermTranslator>(solver_);

  UnorderedTermMap & cache = to_solver_->get_cache();
  Term ns;
  for (auto s : ts_.statevars()) {
    // common variables are next states, unless used for refinement in IC3IA
    // then will refer to current state variables after untiming
    // need to cache both
    cache[to_interpolator_->transfer_term(s)] = s;
    ns = next(s);
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

  for (auto s : free_symbols) {
    assert(s->is_symbol());
    if (s->is_symbolic_const()) {
      // ignore constants
      continue;
    }
    cache[to_interpolator_->transfer_term(s)] = s;
  }
}

ProverResult ModelBasedIC3::check_until(int k)
{
  initialize();
  check_ts();

  // make sure the labels have semantics
  // will check if this has already been done automatically
  set_labels();

  for (int i = 0; i <= k; ++i) {
    ProverResult r = step(i);
    if (r != ProverResult::UNKNOWN) {
      return r;
    }
  }

  return ProverResult::UNKNOWN;
}

bool ModelBasedIC3::witness(vector<UnorderedTermMap> & out)
{
  throw PonoException("ModelBasedIC3 doesn't support witnesses yet.");
}

ProverResult ModelBasedIC3::step(int i)
{
  if (i <= reached_k_) {
    return ProverResult::UNKNOWN;
  }

  if (reached_k_ < 0) {
    return step_0();
  }

  // reached_k_ is the number of transitions that have been checked
  // at this point there are reached_k_ + 1 frames that don't
  // intersect bad, and reached_k_ + 2 frames overall
  assert(reached_k_ + 2 == frames_.size());
  logger.log(1, "Blocking phase at frame {}", i);
  // blocking phase
  while (intersects_bad()) {
    assert(has_proof_goals());
    if (!block_all()) {
      // counter-example
      return ProverResult::FALSE;
    }
  }

  logger.log(1, "Propagation phase at frame {}", i);
  // propagation phase
  push_frame();
  for (size_t j = 1; j < frames_.size() - 1; ++j) {
    if (propagate(j)) {
      assert(j + 1 < frames_.size());
      // save the invariant
      // which is the frame that just had all terms
      // from the previous frames propagated
      set_invar(j + 1);
      return ProverResult::TRUE;
    }
  }

  ++reached_k_;

  return ProverResult::UNKNOWN;
}

ProverResult ModelBasedIC3::step_0()
{
  logger.log(1, "Checking if initial states satisfy property");
  assert(reached_k_ < 0);

  push_solver_context();
  solver_->assert_formula(init_label_);
  solver_->assert_formula(bad_);
  Result r = solver_->check_sat();
  if (r.is_sat()) {
    return ProverResult::FALSE;
  } else {
    assert(r.is_unsat());
    reached_k_ = 0;  // keep reached_k_ aligned with number of frames
  }
  pop_solver_context();
  return ProverResult::UNKNOWN;
}

bool ModelBasedIC3::intersects_bad()
{
  push_solver_context();
  // assert the last frame (conjunction over clauses)
  assert_frame_labels(reached_k_ + 1);
  // see if it intersects with bad
  solver_->assert_formula(bad_);
  Result r = solver_->check_sat();

  if (r.is_sat()) {
    // push bad as a proof goal
    TermVec conjuncts;
    conjunctive_partition(bad_, conjuncts, true);
    Conjunction bad_at_last_frame(solver_, conjuncts);
    add_proof_goal(bad_at_last_frame, reached_k_ + 1, NULL);
  }

  pop_solver_context();

  assert(!r.is_unknown());
  return r.is_sat();
}

bool ModelBasedIC3::get_predecessor(size_t i,
                                    const Conjunction & c,
                                    Conjunction & out_pred)
{
  assert(i > 0);
  assert(i < frames_.size());

  assert(solver_context_ == 0);
  push_solver_context();

  // F[i-1]
  assert_frame_labels(i - 1);
  // -c
  solver_->assert_formula(solver_->make_term(Not, c.term_));
  // Trans
  assert_trans_label();
  // c'
  solver_->assert_formula(next(c.term_));

  Result r = solver_->check_sat();
  if (r.is_sat()) {
    // don't pop now. generalize predecessor needs model values
    // generalize_predecessor will call pop_solver_context();
    out_pred = generalize_predecessor(i, c);
  } else {
    pop_solver_context();

    TermVec assump, red_assump, rem_assump;
    for (auto a : c.conjuncts_) {
      assump.push_back(next(a));
    }

    Term formula = make_and(TermVec{
        get_frame(i - 1), solver_->make_term(Not, c.term_), trans_label_ });

    // filter using unsatcore
    reduce_assump_unsatcore(formula, assump, red_assump, &rem_assump);
    // get current version of red_assump
    TermVec cur_red_assump, cur_rem_assump;
    for (auto a : red_assump) {
      cur_red_assump.push_back(ts_.curr(a));
    }
    for (auto a : rem_assump) {
      cur_rem_assump.push_back(ts_.curr(a));
    }

    fix_if_intersects_initial(cur_red_assump, cur_rem_assump);
    assert(cur_red_assump.size() > 0);
    out_pred = Conjunction(solver_, cur_red_assump);
  }

  assert(!r.is_unknown());
  return r.is_sat();
}

ProofGoal ModelBasedIC3::get_next_proof_goal()
{
  assert(has_proof_goals());
  ProofGoal pg = proof_goals_.back();
  proof_goals_.pop_back();
  return pg;
}

void ModelBasedIC3::add_proof_goal(const Conjunction & c,
                                   size_t i,
                                   shared_ptr<ProofGoal> n)
{
  proof_goals_.push_back(ProofGoal(c, i, n));
}

bool ModelBasedIC3::block_all()
{
  while (has_proof_goals()) {
    ProofGoal pg = get_next_proof_goal();
    // block can fail, which just means a
    // new proof goal will be added
    if (!block(pg) && !pg.idx) {
      // if a proof goal cannot be blocked at zero
      // then there's a counterexample
      return false;
    }
  }
  assert(!has_proof_goals());
  return true;
}

bool ModelBasedIC3::block(const ProofGoal & pg)
{
  const Conjunction & c = pg.conj;
  size_t i = pg.idx;

  logger.log(
      3, "Attempting to block proof goal <{}, {}>", c.term_->to_string(), i);

  assert(i < frames_.size());
  assert(i >= 0);
  // TODO: assert c -> frames_[i]

  if (i == 0) {
    // can't block anymore -- this is a counterexample
    return false;
  }

  Conjunction pred;  // populated by get_predecessor
  if (!get_predecessor(i, c, pred)) {
    // can block this cube
    Term gen_blocking_term = inductive_generalization(i, pred);
    // pred is a subset of c
    logger.log(3, "Blocking term at frame {}: {}", i, c.term_->to_string());
    logger.log(3, " with {}", gen_blocking_term->to_string());

    // if using interpolants, can't count on blocking term being a clause
    // simple heuristic is to get conjunctive partition
    TermVec conjuncts;
    conjunctive_partition(gen_blocking_term, conjuncts, true);
    size_t min_idx = frames_.size();
    for (auto bt : conjuncts) {
      // TODO: fix name -- might not be a clause anymore
      // try to push
      size_t idx = push_blocking_clause(i, bt);
      constrain_frame(idx, bt);
      if (idx < min_idx) {
        min_idx = idx;
      }
    }

    // we're limited by the minimum index that a conjunct could be pushed to
    if (min_idx + 1 < frames_.size()) {
      add_proof_goal(c, min_idx + 1, pg.next);
    }
    return true;
  } else {
    add_proof_goal(pred, i - 1, make_shared<ProofGoal>(pg));
    return false;
  }
}

bool ModelBasedIC3::propagate(size_t i)
{
  assert(i + 1 < frames_.size());

  unordered_set<size_t> indices_to_remove;
  const TermVec & Fi = frames_.at(i);

  push_solver_context();
  assert_frame_labels(i);
  assert_trans_label();

  for (size_t j = 0; j < Fi.size(); ++j) {
    const Term & t = Fi.at(j);

    // Relative inductiveness check
    // Check F[i] /\ t /\ T /\ -t'
    // NOTE: asserting t is redundant because t \in F[i]
    push_solver_context();
    solver_->assert_formula(solver_->make_term(Not, next(t)));

    Result r = solver_->check_sat();
    assert(!r.is_unknown());
    if (r.is_unsat()) {
      // mark for removal
      indices_to_remove.insert(j);
      // add to next frame
      constrain_frame(i + 1, Fi.at(j));
    }

    pop_solver_context();
  }

  pop_solver_context();

  // keep invariant that terms are kept at highest frame
  // where they are known to hold
  // need to remove some terms from the current frame
  TermVec new_frame_i;
  new_frame_i.reserve(Fi.size() - indices_to_remove.size());
  for (size_t j = 0; j < Fi.size(); ++j) {
    if (indices_to_remove.find(j) == indices_to_remove.end()) {
      new_frame_i.push_back(Fi.at(j));
    }
  }

  frames_[i] = new_frame_i;

  return new_frame_i.empty();
}

void ModelBasedIC3::push_frame()
{
  assert(frame_labels_.size() == frames_.size());
  // pushes an empty frame
  frame_labels_.push_back(
      solver_->make_symbol("__frame_label_" + std::to_string(frames_.size()),
                           solver_->make_sort(BOOL)));
  frames_.push_back({});
}

void ModelBasedIC3::constrain_frame(size_t i, const Term & constraint)
{
  assert(i < frame_labels_.size());
  assert(frame_labels_.size() == frames_.size());
  solver_->assert_formula(
      solver_->make_term(Implies, frame_labels_.at(i), constraint));
  frames_.at(i).push_back(constraint);
}

Term ModelBasedIC3::inductive_generalization(size_t i, const Conjunction & c)
{
  Term gen_res = solver_->make_term(Not, c.term_);

  if (options_.ic3_indgen_) {
    if (options_.ic3_indgen_mode_ == 0) {
      UnorderedTermSet keep, core_set;
      TermVec bool_assump, tmp, new_tmp, removed, lits;
      split_eq(solver_, c.conjuncts_, lits);

      if (options_.random_seed_ > 0) {
        shuffle(lits.begin(), lits.end(),
                default_random_engine(options_.random_seed_));
      }

      int iter = 0;
      bool progress = true;
      while (iter <= options_.ic3_gen_max_iter_ && lits.size() > 1 &&
             progress) {
        iter = options_.ic3_gen_max_iter_ > 0 ? iter+1 : iter;
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
              solver_->assert_formula(solver_->make_term(Implies, l, next(t)));
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

      gen_res = solver_->make_term(Not, make_and(lits));

    } else if (options_.ic3_indgen_mode_ == 1) {
      TermVec tmp, lits, red_lits;
      for (auto a : c.conjuncts_) {
        tmp.push_back(next(a));
      }
      split_eq(solver_, tmp, lits);

      // ( (frame /\ trans /\ not(c)) \/ init') /\ c' is unsat
      Term formula = make_and(
          { get_frame(i - 1), trans_label_, solver_->make_term(Not, c.term_) });
      formula = solver_->make_term(Or, formula, next(ts_.init()));
      reduce_assump_unsatcore(formula, lits, red_lits);
      gen_res = solver_->make_term(Not, ts_.curr(make_and(red_lits)));

    } else if (options_.ic3_indgen_mode_ == 2) {
      interpolator_->reset_assertions();

      TermVec conjuncts;
      split_eq(solver_, c.conjuncts_, conjuncts);

      // ( (frame /\ trans /\ not(c)) \/ init') /\ c' is unsat
      Term formula = make_and({ get_frame(i - 1),
                                trans_label_,
                                solver_->make_term(Not, make_and(conjuncts)) });
      formula = solver_->make_term(Or, formula, next(ts_.init()));

      Term int_A = to_interpolator_->transfer_term(formula, BOOL);
      // still use c in B
      // only split equalities in A to encourage more general unsat proofs /
      // interpolants
      Term int_B = to_interpolator_->transfer_term(next(c.term_), BOOL);

      Term interp;
      Result r = interpolator_->get_interpolant(int_A, int_B, interp);
      assert(r.is_unsat());

      Term solver_interp = to_solver_->transfer_term(interp);
      gen_res = ts_.curr(solver_interp);

      logger.log(3, "Got interpolant: {}", gen_res);

    } else {
      assert(false);
    }
  }

  assert(only_curr(gen_res));
  assert(!intersects_initial(solver_->make_term(Not, gen_res)));
  return gen_res;
}

Conjunction ModelBasedIC3::generalize_predecessor(size_t i,
                                                  const Conjunction & c)
{
  DisjointSet ds;
  UnorderedTermMap model;
  const UnorderedTermSet & statevars = ts_.statevars();
  TermVec cube_lits;
  cube_lits.reserve(statevars.size());
  TermVec next_lits;
  next_lits.reserve(statevars.size());
  for (auto v : statevars) {
    Term val = solver_->get_value(v);
    cube_lits.push_back(solver_->make_term(Equal, v, val));
    ds.add(v, val);
    assert(ts_.is_curr_var(v));
    assert(model.find(v) == model.end());
    model[v] = val;

    Term nv = next(v);
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
  for (auto v : inputvars) {
    Term val = solver_->get_value(v);
    input_lits.push_back(solver_->make_term(Equal, v, val));
    assert(model.find(v) == model.end());
    model[v] = val;
  }

  // pop the context that was pushed in get_predecessor
  // now that we have the model
  // now safe to make new calls to the solver
  pop_solver_context();

  // if not generalized, the current state assignments in cube_lits
  // are the predecessor
  Conjunction res(solver_, cube_lits);

  assert(i > 0);
  if (i == 1) {
    // don't need to generalize if i == 1
    // the predecessor is an initial state
    assert(intersects_initial(res.term_));
    return res;
  }

  // should never intersect with a frame before F[i-1]
  // otherwise, this predecessor should have been found
  // in a previous step (before a new frame was pushed)
  assert(i > 1);  // so F[i-2] makes sense
  assert(!intersects(res.term_, get_frame(i - 2)));

  if (options_.ic3_pregen_ && !options_.ic3_functional_preimage_) {
    // add congruent equalities to cube_lits
    for (auto v : statevars) {
      Term t = ds.find(v);
      if (t != v) {
        cube_lits.push_back(solver_->make_term(Equal, t, v));
      }
    }

    Term formula = make_and(input_lits);

    if (ts_.is_deterministic()) {
      formula = solver_->make_term(And, formula, trans_label_);
      formula = solver_->make_term(
          And, formula, solver_->make_term(Not, next(c.term_)));
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
      pre_formula = solver_->make_term(And, pre_formula, get_trans());
      pre_formula = solver_->make_term(And, pre_formula,
                                       solver_->make_term(Not, c.term_));
      pre_formula = solver_->make_term(And, pre_formula, next(c.term_));

      formula = solver_->make_term(And, formula,
                                   solver_->make_term(Not, pre_formula));
    }

    TermVec splits, red_cube_lits, rem_cube_lits;
    split_eq(solver_, cube_lits, splits);
    reduce_assump_unsatcore(formula, splits, red_cube_lits, &rem_cube_lits);
    // should need some assumptions
    // formula should not be unsat on its own
    assert(red_cube_lits.size() > 0);

    // update res to the generalization
    res = Conjunction(solver_, red_cube_lits);

  } else if (options_.ic3_pregen_ && options_.ic3_functional_preimage_) {
    assert(ts_.is_deterministic());

    UnorderedTermMap m;
    for (auto v : inputvars) {
      m[v] = model.at(v);
    }
    for (auto v : statevars) {
      Term nv = next(v);
      m[nv] = model.at(nv);
    }

    Term fun_preimage = solver_->substitute(trans_label_, m);
    TermVec conjuncts;
    conjunctive_partition(fun_preimage, conjuncts, true);
    res = Conjunction(solver_, conjuncts);
  }

  // should not intersect a previous frame
  // or else would have found this predecessor earlier
  assert(i > 1);  // so F[i-2] makes sense
  assert(!intersects(res.term_, get_frame(i - 2)));
  return res;
}

bool ModelBasedIC3::intersects(const Term & A, const Term & B)
{
  push_solver_context();
  solver_->assert_formula(A);
  solver_->assert_formula(B);
  Result r = solver_->check_sat();
  pop_solver_context();
  return r.is_sat();
}

bool ModelBasedIC3::intersects_initial(const Term & t)
{
  return intersects(init_label_, t);
}

void ModelBasedIC3::assert_frame_labels(size_t i) const
{
  // never expecting to assert a frame at base context
  assert(solver_context_ > 0);
  assert(frame_labels_.size() == frames_.size());
  Term assump;
  for (size_t j = 0; j < frame_labels_.size(); ++j) {
    assump = frame_labels_[j];
    if (j < i) {
      // optimization: disable the unused constraints
      // by asserting the negated label
      assump = solver_->make_term(Not, assump);
    }
    assert(assump);  // assert that it's non-null
    solver_->assert_formula(assump);
  }
}

Term ModelBasedIC3::get_frame(size_t i) const
{
  if (i == 0) {
    assert(frames_[0].size() == 1);
    return frames_[0][0];
  }

  Term res = solver_true_;
  for (size_t j = i; j < frames_.size(); ++j) {
    for (auto c : frames_[j]) {
      res = solver_->make_term(And, res, c);
    }
  }
  return res;
}

void ModelBasedIC3::assert_trans_label() const
{
  // shouldn't be a scenario where trans is asserted at base context
  assert(solver_context_ > 0);
  solver_->assert_formula(trans_label_);
}

Term ModelBasedIC3::get_trans() const { return ts_.trans(); };

void ModelBasedIC3::fix_if_intersects_initial(TermVec & to_keep,
                                              const TermVec & rem)
{
  if (rem.size() != 0) {
    Term formula = solver_->make_term(And, init_label_, make_and(to_keep));
    reduce_assump_unsatcore(formula, rem, to_keep);
  }
}

size_t ModelBasedIC3::push_blocking_clause(size_t i, Term c)
{
  push_solver_context();
  solver_->assert_formula(c);
  solver_->assert_formula(solver_->make_term(Not, next(c)));
  assert_trans_label();

  Result r;
  size_t j;
  for (j = i; j + 1 < frames_.size(); ++j) {
    push_solver_context();
    assert_frame_labels(j);
    r = solver_->check_sat();
    pop_solver_context();
    if (r.is_sat()) {
      break;
    }
  }

  pop_solver_context();

  return j;
}

void ModelBasedIC3::reduce_assump_unsatcore(const Term & formula,
                                            const TermVec & assump,
                                            TermVec & out_red,
                                            TermVec * out_rem)
{
  TermVec cand_res = assump;
  TermVec bool_assump, tmp_assump;

  // should be at context 0 or else assertions could be polluted
  assert(solver_context_ == 0);

  push_solver_context();
  solver_->assert_formula(formula);

  // exit if the formula is unsat without assumptions.
  Result r = solver_->check_sat();
  if (r.is_unsat()) {
    pop_solver_context();
    return;
  }

  if (options_.random_seed_ > 0) {
    shuffle(cand_res.begin(), cand_res.end(),
            default_random_engine(options_.random_seed_));
  }

  for (auto a : cand_res) {
    Term l = label(a);
    solver_->assert_formula(solver_->make_term(Implies, l, a));
    bool_assump.push_back(l);
  }

  unsigned iter = 0;
  while (iter <= options_.ic3_gen_max_iter_) {
    iter = options_.ic3_gen_max_iter_ > 0 ? iter+1 : iter;
    r = solver_->check_sat_assuming(bool_assump);
    assert(r.is_unsat());

    bool_assump.clear();
    tmp_assump.clear();

    UnorderedTermSet core_set;
    solver_->get_unsat_core(core_set);
    for (auto a : cand_res) {
      Term l = label(a);
      if (core_set.find(l) != core_set.end()) {
        tmp_assump.push_back(a);
        bool_assump.push_back(l);
      } else if (out_rem) {
        out_rem->push_back(a);
      }
    }

    if (tmp_assump.size() == cand_res.size()) {
      break;
    } else {
      cand_res = tmp_assump;
    }
  }

  pop_solver_context();
  // copy the result
  out_red.insert(out_red.end(), cand_res.begin(), cand_res.end());
}

Term ModelBasedIC3::label(const Term & t)
{
  auto it = labels_.find(t);
  if (it != labels_.end()) {
    return labels_.at(t);
  }

  unsigned i = 0;
  Term l;
  while (true) {
    try {
      l = solver_->make_symbol(
          "assump_" + std::to_string(t->hash()) + "_" + std::to_string(i),
          solver_->make_sort(BOOL));
      break;
    }
    catch (IncorrectUsageException & e) {
      ++i;
    }
    catch (SmtException & e) {
      throw e;
    }
  }

  labels_[t] = l;
  return l;
}

Term ModelBasedIC3::make_and(TermVec vec, SmtSolver slv) const
{
  if (!slv) {
    slv = solver_;
  }

  if (vec.size() == 0) {
    return slv->make_term(true);
  }

  // sort the conjuncts
  std::sort(vec.begin(), vec.end(), term_hash_lt);
  Term res = vec[0];
  for (size_t i = 1; i < vec.size(); ++i) {
    res = slv->make_term(And, res, vec[i]);
  }
  return res;
}

Term ModelBasedIC3::make_or(TermVec vec, SmtSolver slv) const
{
  if (!slv) {
    slv = solver_;
  }

  if (vec.size() == 0) {
    return slv->make_term(false);
  }

  // sort the conjuncts
  std::sort(vec.begin(), vec.end(), term_hash_lt);
  Term res = vec[0];
  for (size_t i = 1; i < vec.size(); ++i) {
    res = slv->make_term(Or, res, vec[i]);
  }
  return res;
}

void ModelBasedIC3::push_solver_context()
{
  solver_->push();
  solver_context_++;
}

void ModelBasedIC3::pop_solver_context()
{
  assert(solver_context_ > 0);
  solver_->pop();
  solver_context_--;
}

void ModelBasedIC3::set_labels()
{
  assert(solver_context_ == 0);  // expecting to be at base context level
  // set semantics of labels
  Sort boolsort = solver_->make_sort(BOOL);
  if (!init_label_) {
    // frame 0 label is identical to init label
    init_label_ = frame_labels_[0];
  }
  if (!trans_label_) {
    trans_label_ = solver_->make_symbol("__trans_label", boolsort);
    solver_->assert_formula(solver_->make_term(Implies, trans_label_, ts_.trans()));
  }
}

void ModelBasedIC3::check_ts() const
{
  // check if there are arrays or uninterpreted sorts and fail if so
  for (auto vec : { ts_.statevars(), ts_.inputvars() }) {
    for (auto st : vec) {
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

bool ModelBasedIC3::only_curr(Term & t) { return ts_.only_curr(t); }

Term ModelBasedIC3::next(const Term & t) const { return ts_.next(t); }

void ModelBasedIC3::set_invar(size_t i) { invar_ = get_frame(i); }

DisjointSet::DisjointSet() {}

DisjointSet::~DisjointSet() {}

void DisjointSet::add(const Term & a, const Term & b)
{
  if (leader_.find(a) != leader_.end()) {
    Term leadera = leader_.at(a);
    UnorderedTermSet & groupa = group_.at(leadera);

    if (leader_.find(b) != leader_.end()) {
      Term leaderb = leader_.at(b);

      if (leadera != leaderb) {
        UnorderedTermSet & groupb = group_.at(leaderb);

        if (leadera <= leaderb) {
          groupa.insert(groupb.begin(), groupb.end());

          for (const Term & t : groupb) {
            leader_[t] = leadera;
          }
          groupb.clear();
          group_.erase(leaderb);

        } else {
          groupb.insert(groupa.begin(), groupa.end());

          for (const Term & t : groupa) {
            leader_[t] = leaderb;
          }
          groupa.clear();
          group_.erase(leadera);
        }
      }

    } else {
      groupa.insert(b);
      leader_[b] = leadera;
    }

  } else if (leader_.find(b) != leader_.end()) {
    Term leaderb = leader_.at(b);
    group_[leaderb].insert(a);
    leader_[a] = leaderb;

  } else {
    if (!a->is_value()) {
      leader_[a] = a;
      leader_[b] = a;
      group_[a] = UnorderedTermSet({ a, b });
    } else {
      assert(!b->is_value());
      leader_[a] = b;
      leader_[b] = b;
      group_[b] = UnorderedTermSet({ a, b });
    }
  }
}

Term DisjointSet::find(const Term & t)
{
  assert(leader_.find(t) != leader_.end());
  return leader_.at(t);
}

void DisjointSet::clear()
{
  leader_.clear();
  group_.clear();
}

}  // namespace pono
