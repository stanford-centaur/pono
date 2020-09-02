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
** \brief Simple implementation of IC3 operating on a functional
**        transition system (exploiting this structure for
**        predecessor computation) and uses models.
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
      true_(solver_->make_term(true)),
      false_(solver_->make_term(false))
{
  // super sets other options
  solver_->set_opt("produce-unsat-cores", "true");
  initialize();
}

ModelBasedIC3::ModelBasedIC3(Property & p, const SmtSolver & slv)
    : super(p, slv),
      true_(solver_->make_term(true)),
      false_(solver_->make_term(false))
{
  initialize();
}

ModelBasedIC3::ModelBasedIC3(const PonoOptions & opt,
                             Property & p,
                             const SolverEnum se)
    : super(opt, p, se),
      true_(solver_->make_term(true)),
      false_(solver_->make_term(false))
{
  // super sets other options
  solver_->set_opt("produce-unsat-cores", "true");
  initialize();
}

ModelBasedIC3::ModelBasedIC3(const PonoOptions & opt,
                             Property & p,
                             const SmtSolver & slv)
    : super(opt, p, slv),
      true_(solver_->make_term(true)),
      false_(solver_->make_term(false))
{
  initialize();
}

ModelBasedIC3::~ModelBasedIC3() {}

void ModelBasedIC3::initialize()
{
  super::initialize();
  frames_.clear();
  proof_goals_.clear();
  // first frame is always the initial states
  frames_.push_back({ ts_.init() });
  push_frame();

  // check if there are arrays and fail if so
  for (auto vec : { ts_.statevars(), ts_.inputvars() }) {
    for (auto st : vec) {
      if (st->get_sort()->get_sort_kind() == ARRAY) {
        throw PonoException("ModelBasedIC3 does not support arrays yet");
      }
    }
  }

  // find UF applications and keep them so they can be
  // included in models
  TermOpCollector toc(solver_);
  UnorderedTermSet uf_apps;
  toc.find_matching_terms(ts_.init(), { Apply }, uf_apps);
  toc.find_matching_terms(ts_.trans(), { Apply }, uf_apps);
  toc.find_matching_terms(bad_, { Apply }, uf_apps);
  extra_model_terms_.clear();
  extra_model_terms_.reserve(uf_apps.size());
  extra_model_terms_.insert(
      extra_model_terms_.begin(), uf_apps.begin(), uf_apps.end());

  // save the uninterpreted functions themselves
  UnorderedTermSet free_symbols;
  get_free_symbols(bad_, free_symbols);
  get_free_symbols(ts_.init(), free_symbols);
  get_free_symbols(ts_.trans(), free_symbols);
  for (auto s : free_symbols) {
    if (s->get_sort()->get_sort_kind() == FUNCTION) {
      ufs_.insert(s);
    }
  }
}

ProverResult ModelBasedIC3::check_until(int k)
{
  for (int i = 0; i <= k; ++i) {
    ProverResult r = step(i);
    if (r != ProverResult::UNKNOWN) {
      return r;
    }
  }

  return ProverResult::UNKNOWN;
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

  solver_->push();
  solver_->assert_formula(ts_.init());
  solver_->assert_formula(bad_);
  Result r = solver_->check_sat();
  if (r.is_sat()) {
    return ProverResult::FALSE;
  } else {
    assert(r.is_unsat());
    reached_k_ = 0;  // keep reached_k_ aligned with number of frames
  }
  solver_->pop();
  return ProverResult::UNKNOWN;
}

bool ModelBasedIC3::intersects_bad()
{
  solver_->push();
  // assert the last frame (conjunction over clauses)
  assert_frame(reached_k_ + 1);
  // see if it intersects with bad
  solver_->assert_formula(bad_);
  Result r = solver_->check_sat();

  if (r.is_sat()) {
    // push bad as a proof goal
    TermVec conjuncts;
    conjunctive_partition(bad_, conjuncts);
    Conjunction bad_at_last_frame(solver_, conjuncts);
    add_proof_goal(bad_at_last_frame, reached_k_ + 1);
  }

  solver_->pop();

  assert(!r.is_unknown());
  return r.is_sat();
}

bool ModelBasedIC3::get_predecessor(size_t i,
                                    const Conjunction & c,
                                    Conjunction & out_pred)
{
  assert(i > 0);
  assert(i < frames_.size());

  solver_->push();

  // F[i-1]
  assert_frame(i - 1);
  // -c
  solver_->assert_formula(solver_->make_term(Not, c.term_));
  // Trans
  solver_->assert_formula(ts_.trans());
  // c'
  solver_->assert_formula(ts_.next(c.term_));

  Result r = solver_->check_sat();
  if (r.is_sat()) {
    // don't pop now. generalize predecessor needs model values
    // generalize_predecessor will call solver_->pop();
    out_pred = generalize_predecessor(i, c);
  } else {
    solver_->pop();

    TermVec assump, red_assump, rem_assump;
    for (auto a : c.conjuncts_) {
      assump.push_back(ts_.next(a));
    }

    Term formula = make_and(TermVec{get_frame(i - 1),
                                    solver_->make_term(Not, c.term_),
                                    ts_.trans()});

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

void ModelBasedIC3::add_proof_goal(const Conjunction & c, size_t i)
{
  proof_goals_.push_back(ProofGoal(c, i));
}

bool ModelBasedIC3::block_all()
{
  while (has_proof_goals()) {
    ProofGoal pg = get_next_proof_goal();
    // block can fail, which just means a
    // new proof goal will be added
    if (!block(pg) && !pg.second) {
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
  const Conjunction & c = pg.first;
  size_t i = pg.second;

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
    conjunctive_partition(gen_blocking_term, conjuncts);
    size_t min_idx = frames_.size();
    for (auto bt : conjuncts) {
      // TODO: fix name -- might not be a clause anymore
      // try to push
      size_t idx = push_blocking_clause(i, bt);
      frames_[idx].push_back(bt);
      if (idx < min_idx) {
        min_idx = idx;
      }
    }

    // we're limited by the minimum index that a conjunct could be pushed to
    if (min_idx + 1 < frames_.size()) {
      add_proof_goal(c, min_idx + 1);
    }
    return true;
  } else {
    add_proof_goal(pred, i - 1);
    return false;
  }
}

bool ModelBasedIC3::propagate(size_t i)
{
  assert(i + 1 < frames_.size());

  unordered_set<size_t> indices_to_remove;
  const TermVec & Fi = frames_.at(i);

  solver_->push();
  assert_frame(i);
  solver_->assert_formula(ts_.trans());

  for (size_t j = 0; j < Fi.size(); ++j) {
    const Term & t = Fi[j];

    // Relative inductiveness check
    // Check F[i] /\ t /\ T /\ -t'
    // NOTE: asserting t is redundant because t \in F[i]
    solver_->push();
    solver_->assert_formula(solver_->make_term(Not, ts_.next(t)));

    Result r = solver_->check_sat();
    assert(!r.is_unknown());
    if (r.is_unsat()) {
      // mark for removal
      indices_to_remove.insert(j);
      // add to next frame
      frames_[i + 1].push_back(Fi[j]);
    }

    solver_->pop();
  }

  solver_->pop();

  // keep invariant that terms are kept at highest frame
  // where they are known to hold
  // need to remove some terms from the current frame
  TermVec new_frame_i;
  new_frame_i.reserve(Fi.size() - indices_to_remove.size());
  for (size_t j = 0; j < Fi.size(); ++j) {
    if (indices_to_remove.find(j) == indices_to_remove.end()) {
      new_frame_i.push_back(Fi[j]);
    }
  }

  frames_[i] = new_frame_i;

  return new_frame_i.empty();
}

void ModelBasedIC3::push_frame()
{
  // pushes an empty frame
  frames_.push_back({});
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
            solver_->push();
            assert_frame(i - 1);
            solver_->assert_formula(ts_.trans());
            solver_->assert_formula(solver_->make_term(Not, tmp_and_term));

            Term l;
            bool_assump.clear();
            for (auto t : tmp) {
              l = label(t);
              solver_->assert_formula(solver_->make_term(Implies,
                                                         l, ts_.next(t)));
              bool_assump.push_back(l);
            }

            Result r = solver_->check_sat_assuming(bool_assump);
            assert(!r.is_unknown());

            if (r.is_sat()) {
              // we cannot drop a
              solver_->pop();
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

              solver_->pop();

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
        tmp.push_back(ts_.next(a));
      }
      split_eq(solver_, tmp, lits);

      // ( (frame /\ trans /\ not(c)) \/ init') /\ c' is unsat
      Term formula = make_and({get_frame(i - 1), ts_.trans(),
                               solver_->make_term(Not, c.term_)});
      formula = solver_->make_term(Or, formula, ts_.next(ts_.init()));
      reduce_assump_unsatcore(formula, lits, red_lits);
      gen_res = solver_->make_term(Not, ts_.curr(make_and(red_lits)));

    } else if (options_.ic3_indgen_mode_ == 2) {
      // ( (frame /\ trans /\ not(c)) \/ init') /\ c' is unsat
      Term formula = make_and({get_frame(i - 1), ts_.trans(),
                               solver_->make_term(Not, c.term_)});
      formula = solver_->make_term(Or, formula, ts_.next(ts_.init()));

      // TODO make an option to set the interpolator
      SmtSolver interpolator = create_interpolating_solver(MSAT_INTERPOLATOR);
      TermTranslator to_interpolator(interpolator);
      TermTranslator to_solver(solver_);

      UnorderedTermMap & cache = to_solver.get_cache();
      Term ns;
      for (auto s : ts_.statevars()) {
        // common variables will be next state
        // so that's all we need
        // better not to cache the others, now if there's a bug
        // where the shared variables are not respected, the term
        // translator will throw an exception
        ns = ts_.next(s);
        cache[to_interpolator.transfer_term(ns)] = ns;
      }

      // need to copy over UF as well
      for (auto uf : ufs_) {
        assert(uf->get_sort()->get_sort_kind() == FUNCTION);
        cache[to_interpolator.transfer_term(uf)] = uf;
      }

      Term int_A = to_interpolator.transfer_term(formula, BOOL);
      Term int_B = to_interpolator.transfer_term(ts_.next(c.term_), BOOL);

      Term interp;
      Result r = interpolator->get_interpolant(int_A, int_B, interp);
      assert(r.is_unsat());

      Term solver_interp = to_solver.transfer_term(interp);
      gen_res = ts_.curr(solver_interp);
      assert(ts_.only_curr(gen_res));

      logger.log(3, "Got interpolant: {}", gen_res);

    } else {
      assert(false);
    }
  }

  assert(!intersects_initial(solver_->make_term(Not, gen_res)));
  return gen_res;
}

Conjunction ModelBasedIC3::generalize_predecessor(size_t i,
                                                  const Conjunction & c)
{
  const UnorderedTermSet & statevars = ts_.statevars();
  TermVec cube_lits;
  DisjointSet ds;
  cube_lits.reserve(statevars.size() + extra_model_terms_.size());
  for (auto v : statevars) {
    Term val = solver_->get_value(v);
    cube_lits.push_back(solver_->make_term(Equal, v, val));
    ds.add(v, val);
  }

  // collect input assignments
  UnorderedTermMap input_assignments;
  const UnorderedTermSet & inputvars = ts_.inputvars();
  TermVec input_lits;
  input_lits.reserve(inputvars.size());
  for (auto v : inputvars) {
    Term val = solver_->get_value(v);
    input_assignments[v] = val;
    input_lits.push_back(solver_->make_term(Equal, v, val));
  }

  // get other important model values
  for (auto t : extra_model_terms_) {
    Term val = solver_->get_value(t);
    Term t_subs = solver_->substitute(t, input_assignments);
    if (ts_.only_curr(t_subs)) {
      cube_lits.push_back(solver_->make_term(Equal, t_subs, val));
      ds.add(t, val);
    }
  }

  Conjunction res(solver_, cube_lits);

  if (options_.ic3_cexgen_ && !options_.ic3_functional_preimage_) {
    // add congruent equalities to cube_lits
    for (auto v : statevars) {
      Term t = ds.find(v);
      if (t != v) {
        cube_lits.push_back(solver_->make_term(Equal, t, v));
      }
    }

    // TODO: figure out if anything else should be done for UF

    // collect next statevars assignments
    TermVec next_lits;
    if (!ts_.is_deterministic()) {
      next_lits.reserve(statevars.size());
      for (auto v : statevars) {
        Term nv = ts_.next(v);
        next_lits.push_back(solver_->make_term(Equal, nv,
                                               solver_->get_value(nv)));
      }

      for (auto t : extra_model_terms_) {
        if (ts_.only_curr(t)) {
          Term nt = ts_.next(t);
          next_lits.push_back(
              solver_->make_term(Equal, nt, solver_->get_value(nt)));
        } else {
          solver_->make_term(Equal, t, solver_->get_value(t));
        }
      }
    }

    solver_->pop();

    Term formula = make_and(input_lits);

    if (ts_.is_deterministic()) {
      formula = solver_->make_term(And, formula, ts_.trans());
      formula = solver_->make_term(And, formula,
                                   solver_->make_term(Not, ts_.next(c.term_)));
    } else {
      formula = solver_->make_term(And, formula, make_and(next_lits));

      // preimage formula
      Term pre_formula = get_frame(i - 1);
      pre_formula = solver_->make_term(And, pre_formula, ts_.trans());
      pre_formula = solver_->make_term(And, pre_formula,
                                       solver_->make_term(Not, c.term_));
      pre_formula = solver_->make_term(And, pre_formula, ts_.next(c.term_));

      formula = solver_->make_term(And, formula,
                                   solver_->make_term(Not, pre_formula));
    }

    TermVec splits, red_cube_lits;
    split_eq(solver_, cube_lits, splits);
    reduce_assump_unsatcore(formula, splits, red_cube_lits);
    assert(red_cube_lits.size() > 0);
    res = Conjunction(solver_, red_cube_lits);

  } else if (options_.ic3_cexgen_ && options_.ic3_functional_preimage_) {
    assert(ts_.is_functional());

    UnorderedTermMap m;
    const UnorderedTermSet & inputvars = ts_.inputvars();
    for (auto v : inputvars) {
      m[v] = solver_->get_value(v);
    }
    for (auto v : statevars) {
      Term nv = ts_.next(v);
      m[nv] = solver_->get_value(nv);
    }

    solver_->pop();

    Term fun_preimage = solver_->substitute(ts_.trans(), m);
    TermVec conjuncts;
    conjunctive_partition(fun_preimage, conjuncts);
    res = Conjunction(solver_, conjuncts);

  } else {
    solver_->pop();
  }

  return res;
}

bool ModelBasedIC3::intersects_initial(const Term & t) const
{
  solver_->push();
  solver_->assert_formula(t);
  solver_->assert_formula(ts_.init());
  Result r = solver_->check_sat();
  solver_->pop();
  return r.is_sat();
}

void ModelBasedIC3::assert_frame(size_t i) const
{
  solver_->assert_formula(get_frame(i));
}

Term ModelBasedIC3::get_frame(size_t i) const
{
  if (i == 0) {
    assert(frames_[0].size() == 1);
    return frames_[0][0];
  }

  Term res = true_;
  for (size_t j = i; j < frames_.size(); ++j) {
    for (auto c : frames_[j]) {
      res = solver_->make_term(And, res, c);
    }
  }
  return res;
}

void ModelBasedIC3::fix_if_intersects_initial(TermVec & to_keep,
                                              const TermVec & rem)
{
  if (rem.size() != 0) {
    Term formula = solver_->make_term(And, ts_.init(), make_and(to_keep));
    reduce_assump_unsatcore(formula, rem, to_keep);
  }
}

size_t ModelBasedIC3::push_blocking_clause(size_t i, Term c)
{
  solver_->push();
  solver_->assert_formula(c);
  solver_->assert_formula(solver_->make_term(Not, ts_.next(c)));
  solver_->assert_formula(ts_.trans());

  Result r;
  size_t j;
  for (j = i; j + 1 < frames_.size(); ++j) {
    solver_->push();
    assert_frame(j);
    r = solver_->check_sat();
    solver_->pop();
    if (r.is_sat()) {
      break;
    }
  }

  solver_->pop();

  return j;
}

void ModelBasedIC3::reduce_assump_unsatcore(const Term & formula,
                                            const TermVec & assump,
                                            TermVec & out_red,
                                            TermVec * out_rem)
{
  TermVec cand_res = assump;
  TermVec bool_assump, tmp_assump;

  solver_->push();
  solver_->assert_formula(formula);

  // exit if the formula is unsat without assumptions.
  Result r = solver_->check_sat();
  if (r.is_unsat()) {
    solver_->pop();
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

  solver_->pop();
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

Term ModelBasedIC3::make_and(smt::TermVec vec) const
{
  if (vec.size() == 0) {
    return true_;
  }

  // sort the conjuncts
  std::sort(vec.begin(), vec.end(), term_hash_lt);
  Term res = vec[0];
  for (size_t i = 1; i < vec.size(); ++i) {
    res = solver_->make_term(And, res, vec[i]);
  }
  return res;
}

Term ModelBasedIC3::make_or(smt::TermVec vec) const
{
  if (vec.size() == 0) {
    return false_;
  }

  // sort the conjuncts
  std::sort(vec.begin(), vec.end(), term_hash_lt);
  Term res = vec[0];
  for (size_t i = 1; i < vec.size(); ++i) {
    res = solver_->make_term(Or, res, vec[i]);
  }
  return res;
}

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

}  // namespace pono
