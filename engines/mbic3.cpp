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
#include "engines/mbic3.h"

#include "utils/logger.h"

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
      if (s->get_sort_kind() == SortKind::BOOL) {
        out.push_back(a);  // assert(false);
      } else if (s->get_sort_kind() == SortKind::BV) {
        out.push_back(solver->make_term(BVUle, args[0], args[1]));
        out.push_back(solver->make_term(BVUle, args[1], args[0]));
      } else {
        out.push_back(solver->make_term(Le, args[0], args[1]));
        out.push_back(solver->make_term(Le, args[1], args[0]));
      }
    } else {
      out.push_back(a);
    }
  }
}

// Clause / Cube implementations

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

Clause::Clause(const smt::SmtSolver & solver, const smt::TermVec & lits)
    : lits_(lits)
{
  // sort literals
  std::sort(lits_.begin(), lits_.end(), term_hash_lt);
  // shouldn't have an empty clause
  assert(lits_.size());

  // create term
  term_ = lits_[0];
  for (size_t i = 1; i < lits_.size(); ++i) {
    term_ = solver->make_term(smt::Or, term_, lits_[i]);
  }
}

Cube::Cube(const smt::SmtSolver & solver, const smt::TermVec & lits)
    : lits_(lits)
{
  // sort literals
  std::sort(lits_.begin(), lits_.end(), term_hash_lt);
  // shouldn't have an empty cube
  assert(lits_.size());

  // create term
  term_ = lits_[0];
  for (size_t i = 1; i < lits_.size(); ++i) {
    term_ = solver->make_term(smt::And, term_, lits_[i]);
  }
}

// main implementations

ModelBasedIC3::ModelBasedIC3(const Property & p, SolverEnum se)
    : super(p, se),
      true_(solver_->make_term(true)),
      false_(solver_->make_term(false))
{
  initialize();
}

ModelBasedIC3::ModelBasedIC3(const Property & p, const SmtSolver & slv)
    : super(p, slv),
      true_(solver_->make_term(true)),
      false_(solver_->make_term(false))
{
  initialize();
}

ModelBasedIC3::ModelBasedIC3(const PonoOptions & opt,
                             const Property & p,
                             const SolverEnum se)
    : super(opt, p, se),
      true_(solver_->make_term(true)),
      false_(solver_->make_term(false))
{
  initialize();
}

ModelBasedIC3::ModelBasedIC3(const PonoOptions & opt,
                             const Property & p,
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
    assert(!proof_goals_.empty());
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
    // create a proof goal for the bad state
    const UnorderedTermSet & statevars = ts_.statevars();
    TermVec cube_vec;
    cube_vec.reserve(statevars.size());
    Term eq;
    for (auto sv : statevars) {
      eq = solver_->make_term(Equal, sv, solver_->get_value(sv));
      cube_vec.push_back(eq);
    }

    // reduce the bad cube
    solver_->pop();
    solver_->push();

    Term neg_assert =
      solver_->make_term(Not, solver_->make_term(And, bad_,
                                                 get_frame(reached_k_ + 1)));
    solver_->assert_formula(neg_assert);

    TermVec splits;
    split_eq(solver_, cube_vec, splits);

    TermVec bool_assump;
    for (auto a : splits) {
      unsigned i = 0;
      Term b = label(a);
      bool_assump.push_back(b);
      solver_->assert_formula(solver_->make_term(Implies, b, a));
    }

    Result rr = solver_->check_sat_assuming(bool_assump);
    assert(rr.is_unsat());
    TermVec core = solver_->get_unsat_core();
    UnorderedTermSet core_set(core.begin(), core.end());
    TermVec red_lits;
    for (size_t j = 0; j < bool_assump.size(); ++j) {
      if (core_set.find(bool_assump[j]) != core_set.end()) {
        red_lits.push_back(splits[j]);
      }
    }

    Cube c(solver_, red_lits);
    proof_goals_[reached_k_ + 1].push_back(c);
  }

  solver_->pop();

  assert(!r.is_unknown());
  return r.is_sat();
}

bool ModelBasedIC3::get_predecessor(size_t i, const Cube & c, Cube & out_pred)
{
  solver_->push();
  assert(i > 0);
  assert(i < frames_.size());
  // F[i-1]
  assert_frame(i - 1);
  // -c
  solver_->assert_formula(solver_->make_term(Not, c.term_));
  // Trans
  solver_->assert_formula(ts_.trans());

  // c'
  TermVec bool_assump;
  for (auto a : c.lits_) {
    unsigned i = 0;
    Term b = label(a);
    bool_assump.push_back(b);
    solver_->assert_formula(solver_->make_term(Implies, b, ts_.next(a)));
  }

  Result r = solver_->check_sat_assuming(bool_assump);
  if (r.is_sat()) {
    // don't pop now. generalize predecessor needs model values
    // generalize_predecessor will call solver_->pop();
    out_pred = generalize_predecessor(i, c);
  } else {
    // filter using unsatcore
    TermVec red_lits, rem;
    TermVec core = solver_->get_unsat_core();
    UnorderedTermSet core_set(core.begin(), core.end());
    for (size_t j = 0; j < bool_assump.size(); ++j) {
      if (core_set.find(bool_assump[j]) != core_set.end()) {
        red_lits.push_back(c.lits_[j]);
      } else {
        rem.push_back(c.lits_[j]);
      }
    }

    solver_->pop();

    fix_if_intersects_initial(red_lits, rem);
    out_pred = Cube(solver_, red_lits);
  }

  assert(!r.is_unknown());
  return r.is_sat();
}

ProofGoal ModelBasedIC3::get_next_proof_goal()
{
  assert(has_proof_goals());
  size_t i = proof_goals_.begin()->first;
  ProofGoal pg(proof_goals_.at(i).back(), i);
  // need to remove the proof goal and the map
  // entry if there are no more cubes in the vector
  proof_goals_.at(i).pop_back();
  if (proof_goals_.at(i).empty()) {
    proof_goals_.erase(i);
  }
  return pg;
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
  const Cube & c = pg.first;
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

  Cube pred;  // populated by get_predecessor
  if (!get_predecessor(i, c, pred)) {
    // can block this cube
    Term gen_blocking_term = inductive_generalization(i, pred);
    //pred is a subset of c
    logger.log(3, "Blocking term at frame {}: {}", i, c.term_->to_string());
    logger.log(3, " with {}", gen_blocking_term->to_string());
    frames_[i].push_back(gen_blocking_term);
    return true;
  } else {
    // add a new proof goal
    proof_goals_[i - 1].push_back(pred);
    return false;
  }
}

bool ModelBasedIC3::propagate(size_t i)
{
  assert(i + 1 < frames_.size());

  unordered_set<size_t> indices_to_remove;
  TermVec & Fi = frames_.at(i);

  solver_->push();
  assert_frame(i);
  solver_->assert_formula(ts_.trans());

  for (size_t j = 0; j < Fi.size(); ++j) {
    Term & t = Fi[j];

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

Term ModelBasedIC3::inductive_generalization(size_t i, const Cube & c)
{
  Term res_cube = make_and(c.lits_);

  if (options_.ic3_indgen_) {
    bool progress = true;
    TermVec bool_assump, tmp, new_tmp;
    UnorderedTermSet keep;
    TermVec lits;
    split_eq(solver_, c.lits_, lits);

    int iter = 0;
    // max 2 iterations
    while (++iter <= 2 && lits.size() > 1 && progress) {
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

          bool_assump.clear();
          for (auto t : tmp) {
            Term l = label(t);
            solver_->assert_formula(
                solver_->make_term(Implies, l, ts_.next(t)));
            bool_assump.push_back(l);
          }

          Result r = solver_->check_sat_assuming(bool_assump);
          assert(!r.is_unknown());

          if (r.is_sat()) {
            // we cannot drop a
            solver_->pop();
          } else {
            // filter using unsatcore
            TermVec removed;
            new_tmp.clear();
            TermVec core = solver_->get_unsat_core();
            UnorderedTermSet core_set(core.begin(), core.end());
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
            break; // next iteration
          }
        }
      }

      progress = lits.size() < prev_size;
    }

    res_cube = make_and(lits);
  }

  assert(!intersects_initial(res_cube));
  return solver_->make_term(Not, res_cube);
}

Clause ModelBasedIC3::down(size_t i, const Clause & c) const
{
  // TODO: implement this when implementing inductive generalization
  // For now, just a stub
  throw PonoException("Not yet implemented");
}

Cube ModelBasedIC3::generalize_predecessor(size_t i, const Cube & c)
{
  const UnorderedTermSet & statevars = ts_.statevars();
  TermVec cube_lits;
  cube_lits.reserve(statevars.size());
  for (auto v : statevars) {
    cube_lits.push_back(solver_->make_term(Equal, v, solver_->get_value(v)));
  }

  Cube res = Cube(solver_, cube_lits);

  if (ts_.is_functional() && options_.ic3_cexgen_) {
    // collect input assignments
    const UnorderedTermSet & inputvars = ts_.inputvars();
    TermVec input_lits;
    input_lits.reserve(inputvars.size());
    for (auto v : inputvars) {
      input_lits.push_back(solver_->make_term(Equal, v, solver_->get_value(v)));
    }

    solver_->pop();
    solver_->push();
    // assert input assignments
    for (auto & a : input_lits) {
      solver_->assert_formula(a);
    }
    // assert T
    solver_->assert_formula(ts_.trans());
    // assert -c'
    solver_->assert_formula(solver_->make_term(Not, ts_.next(c.term_)));

    // make and assert assumptions
    TermVec bool_assump;
    for (auto a : cube_lits) {
      unsigned i = 0;
      Term b = label(a);
      bool_assump.push_back(b);
      solver_->assert_formula(solver_->make_term(Implies, b, a));
    }

    Result r = solver_->check_sat_assuming(bool_assump);
    assert(r.is_unsat());

    // filter using unsatcore
    TermVec red_cube_lits;
    TermVec core = solver_->get_unsat_core();
    UnorderedTermSet core_set(core.begin(), core.end());
    for (size_t j = 0; j < bool_assump.size(); ++j) {
      if (core_set.find(bool_assump[j]) != core_set.end()) {
        red_cube_lits.push_back(cube_lits[j]);
      }
    }

    assert(red_cube_lits.size() > 0);
    res = Cube(solver_, red_cube_lits);
  } else {
    // TODO: write generalize_pred for relational transition systemms
  }

  solver_->pop();

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
  if (rem.size() == 0) {
    return;
  }

  solver_->push();
  solver_->assert_formula(ts_.init());
  for (auto a : to_keep) {
    solver_->assert_formula(a);
  }

  Result r = solver_->check_sat();
  if (r.is_sat()) {
    TermVec bool_assump;
    for (auto a : rem) {
      Term l = label(a);
      solver_->assert_formula(solver_->make_term(Implies, l, a));
      bool_assump.push_back(l);
    }

    r = solver_->check_sat_assuming(bool_assump);
    assert(r.is_unsat());

    TermVec core = solver_->get_unsat_core();
    UnorderedTermSet core_set(core.begin(), core.end());
    for (size_t j = 0; j < bool_assump.size(); ++j) {
      if (core_set.find(bool_assump[j]) != core_set.end()) {
        to_keep.push_back(rem[j]);
      }
    }
  }

  solver_->pop();
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

Term ModelBasedIC3::make_and(const smt::TermVec & vec) const
{
  assert(vec.size() > 0);
  Term res = vec[0];
  for (size_t i = 1; i < vec.size(); ++i) {
    res = solver_->make_term(And, res, vec[i]);
  }
  return res;
}

}  // namespace pono
