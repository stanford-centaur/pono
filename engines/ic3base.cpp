/*********************                                                  */
/*! \file ic3base.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan, Florian Lonsing
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Abstract base class implementation of IC3 parameterized by
**        the unit used in frames, pre-image computation, and inductive
**        and predecessor generalization techniques.
**
**/

#include "engines/ic3base.h"

#include <algorithm>

#include "assert.h"
#include "smt/available_solvers.h"
#include "utils/logger.h"

using namespace smt;
using namespace std;

namespace pono {

// helper functions

/** Less than comparison of the hash of two terms
 *  for use in sorting
 *  @param t0 the first term
 *  @param t1 the second term
 *  @return true iff t0's hash is less than t1's hash
 */
static bool term_hash_lt(const smt::Term & t0, const smt::Term & t1)
{
  return (t0->hash() < t1->hash());
}

/** IC3Base */

IC3Base::IC3Base(Property & p, SolverEnum se)
    : super(p, se), reducer_(create_solver(se)), solver_context_(0)
{
}

IC3Base::IC3Base(Property & p, const SmtSolver & s)
    : super(p, s),
      reducer_(create_solver(s->get_solver_enum())),
      solver_context_(0)
{
}

IC3Base::IC3Base(const PonoOptions & opt, Property & p, SolverEnum se)
    : super(opt, p, se), reducer_(create_solver(se)), solver_context_(0)
{
}

IC3Base::IC3Base(const PonoOptions & opt, Property & p, const SmtSolver & s)
    : super(opt, p, s),
      reducer_(create_solver(s->get_solver_enum())),
      solver_context_(0)
{
}

void IC3Base::initialize()
{
  if (initialized_)
  {
    return;
  }

  super::initialize();

  assert(solver_context_ == 0);  // expecting to be at base context level
  solver_true_ = solver_->make_term(true);

  // abstract the transition relation if this is a CEGAR implementation
  // otherwise it is a No-Op
  abstract();

  // ts_ should not only be set but also be the correct one at this point
  // e.g. if abstracting, it now points to the abstract transition system
  assert(ts_);

  // check whether this flavor of IC3 can be applied to this transition system
  check_ts();

  frames_.clear();
  frame_labels_.clear();
  proof_goals_.clear();
  // first frame is always the initial states
  push_frame();
  // can't use constrain_frame for initial states because not guaranteed to be
  // an IC3Formula it's handled specially
  solver_->assert_formula(
      solver_->make_term(Implies, frame_labels_.at(0), ts_->init()));
  push_frame();

  // set semantics of TS labels
  Sort boolsort = solver_->make_sort(BOOL);
  assert(!init_label_);
  assert(!trans_label_);
  // frame 0 label is identical to init label
  init_label_ = frame_labels_[0];
  trans_label_ = solver_->make_symbol("__trans_label", boolsort);
  solver_->assert_formula(
      solver_->make_term(Implies, trans_label_, ts_->trans()));
}

ProverResult IC3Base::check_until(int k)
{
  initialize();
  // make sure derived class implemented initialize and called
  // this version of initialize with super::initialize or
  // (for experts only) set the initialized_ flag without
  // ever initializing base classes
  assert(initialized_);

  ProverResult res;
  RefineResult ref_res;
  int i = reached_k_ + 1;
  assert(i >= 0);
  while (i <= k) {
    // reset cex_pg_ to null
    // there might be multiple abstract traces if there's a derived class
    // doing abstraction refinement
    cex_pg_ = ProofGoal();

    res = step(i);
    ref_res = REFINE_NONE;  // just a default value

    if (res == ProverResult::TRUE) {
      return res;
    } else if (res == ProverResult::FALSE) {
      // expecting cex_pg_ to be non-null and point to the first proof goal in a
      // trace
      assert(cex_pg_.target.term);
      ref_res = refine();
      if (ref_res == RefineResult::REFINE_NONE) {
        // found a concrete counterexample
        return res;
      } else if (ref_res == RefineResult::REFINE_FAIL) {
        logger.log(1, "Failed in refinement.");
        return ProverResult::UNKNOWN;
      }
    }

    // two cases
    // got unknown, so keep going
    // got false, but was able to refine successfully
    assert(res == ProverResult::UNKNOWN
           || (res == ProverResult::FALSE
               && ref_res == RefineResult::REFINE_SUCCESS));

    // increment i, unless there was a refinement step just done
    if (ref_res != RefineResult::REFINE_SUCCESS) {
      i++;
    }
  }

  return ProverResult::UNKNOWN;
}

bool IC3Base::witness(std::vector<smt::UnorderedTermMap> & out)
{
  throw PonoException("IC3 witness NYI");
}

// Protected Methods

IC3Formula IC3Base::ic3formula_disjunction(const TermVec & c) const
{
  assert(c.size());
  Term term = c.at(0);
  for (size_t i = 1; i < c.size(); ++i) {
    term = solver_->make_term(Or, term, c[i]);
  }
  IC3Formula res(term, c, true);
  return res;
}

IC3Formula IC3Base::ic3formula_conjunction(const TermVec & c) const
{
  assert(c.size());
  Term term = c.at(0);
  for (size_t i = 1; i < c.size(); ++i) {
    term = solver_->make_term(And, term, c[i]);
  }
  IC3Formula res(term, c, false);
  return res;
}

IC3Formula IC3Base::ic3formula_negate(const IC3Formula & u) const
{
  const TermVec & children = u.children;
  assert(!u.is_null());
  assert(children.size());

  TermVec neg_children;
  neg_children.reserve(children.size());
  Term nc = smart_not(children.at(0));

  bool is_clause = u.is_disjunction();
  Term term = nc;
  neg_children.push_back(nc);
  for (size_t i = 1; i < children.size(); ++i) {
    nc = smart_not(children[i]);
    neg_children.push_back(nc);
    if (is_clause) {
      // negation is a cube
      term = solver_->make_term(And, term, nc);
    } else {
      // negation is a clause
      term = solver_->make_term(Or, term, nc);
    }
  }
  IC3Formula res(term, neg_children, !is_clause);
  return res;
}

bool IC3Base::intersects_bad()
{
  push_solver_context();
  // assert the last frame (conjunction over clauses)
  assert_frame_labels(reached_k_ + 1);
  // see if it intersects with bad
  solver_->assert_formula(bad_);
  Result r = solver_->check_sat();

  if (r.is_sat()) {
    IC3Formula c = get_model_ic3formula();
    // reduce c
    TermVec red_c;
    reducer_.reduce_assump_unsatcore(smart_not(bad_), c.children, red_c);

    add_proof_goal(ic3formula_conjunction(red_c), reached_k_ + 1, NULL);
  }

  pop_solver_context();

  assert(!r.is_unknown());
  return r.is_sat();
}

ProverResult IC3Base::step(int i)
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
      invar_ = get_frame_term(j + 1);
      return ProverResult::TRUE;
    }
  }

  ++reached_k_;

  return ProverResult::UNKNOWN;
}

ProverResult IC3Base::step_0()
{
  logger.log(1, "Checking if initial states satisfy property");
  assert(reached_k_ < 0);

  push_solver_context();
  solver_->assert_formula(init_label_);
  solver_->assert_formula(bad_);
  Result r = solver_->check_sat();
  if (r.is_sat()) {
    IC3Formula c = get_model_ic3formula();
    cex_pg_ = ProofGoal(c, 0, nullptr);
    pop_solver_context();
    return ProverResult::FALSE;
  } else {
    assert(r.is_unsat());
    reached_k_ = 0;  // keep reached_k_ aligned with number of frames
  }
  pop_solver_context();
  return ProverResult::UNKNOWN;
}

bool IC3Base::rel_ind_check(size_t i,
                            const IC3Formula & c,
                            vector<IC3Formula> & out)
{
  assert(i > 0);
  assert(i < frames_.size());
  // expecting to be the polarity for proof goals, not frames
  // e.g. a conjunction
  assert(!c.is_disjunction());
  assert(!out.size());  // expecting to get an empty vector to populate

  assert(solver_context_ == 0);
  push_solver_context();

  // F[i-1]
  assert_frame_labels(i - 1);
  // -c
  solver_->assert_formula(solver_->make_term(Not, c.term));
  // Trans
  assert_trans_label();
  // c'
  solver_->assert_formula(ts_->next(c.term));

  Result r = solver_->check_sat();
  if (r.is_sat()) {
    IC3Formula predecessor;
    if (options_.ic3_pregen_) {
      predecessor = generalize_predecessor(i, c);
    } else {
      predecessor = get_model_ic3formula();
    }
    assert(ic3formula_check_valid(predecessor));
    out.push_back(predecessor);
    pop_solver_context();
  } else {
    // TODO: consider automatically taking advantage
    //       of an unsat core. Took it out for now (was in MBIC3)
    //       because it needs to work for any IC3Formula
    //       Maybe IC3Formula needs to know how to generalize itself
    //         or at least how to make a conjunctive partition
    //         or it's possible they all can function approximately the same
    //       would also have to move the pop_solver_context later
    pop_solver_context();
    if (options_.ic3_indgen_) {
      assert(solver_context_ == 0); // important that there are no lingering assertions
      out = inductive_generalization(i, c);
    } else {
      IC3Formula blocking_unit = ic3formula_negate(c);
      out.push_back(blocking_unit);
    }
    Term conj = solver_->make_term(true);
    for (auto u : out) {
      conj = solver_->make_term(And, conj, u.term);
      assert(ic3formula_check_valid(u));
      assert(ts_->only_curr(u.term));
    }
    assert(!check_intersects_initial(solver_->make_term(Not, conj)));
  }
  assert(solver_context_ == 0);

  if (r.is_sat()) {
    // for now, assuming that there's only one predecessor produced
    assert(out.size() == 1);
    // this check needs to be here after the solver context has been popped
    // if i == 1 and there's a predecessor, then it should be an initial state
    assert(i != 1 || check_intersects_initial(out.at(0).term));

    // should never intersect with a frame before F[i-1]
    // otherwise, this predecessor should have been found
    // in a previous step (before a new frame was pushed)
    assert(i < 2 || !check_intersects(out.at(0).term, get_frame_term(i - 2)));
  }

  assert(!r.is_unknown());
  return r.is_unsat();
}

// Helper methods

bool IC3Base::block_all()
{
  while (has_proof_goals()) {
    const ProofGoal &pg = get_next_proof_goal();
    if (is_blocked(pg)) {
      logger.log(3, "Skipping already blocked proof goal <{}, {}>",
                 pg.target.term->to_string(), pg.idx);
      continue;
    };
    // block can fail, which just means a
    // new proof goal will be added
    if (!block(pg) && !pg.idx) {
      // if a proof goal cannot be blocked at zero
      // then there's a counterexample
      cex_pg_ = pg;
      return false;
    }
  }
  assert(!has_proof_goals());
  return true;
}

bool IC3Base::block(const ProofGoal & pg)
{
  const IC3Formula & c = pg.target;
  size_t i = pg.idx;

  logger.log(
      3, "Attempting to block proof goal <{}, {}>", c.term->to_string(), i);

  assert(i < frames_.size());
  assert(i >= 0);
  // TODO: assert c -> frames_[i]

  if (i == 0) {
    // can't block anymore -- this is a counterexample
    return false;
  }

  vector<IC3Formula> collateral;  // populated by rel_ind_check
  if (rel_ind_check(i, c, collateral)) {
    // collateral is a vector of blocking units
    assert(collateral.size());
    logger.log(3, "Blocking term at frame {}: {}", i, c.term->to_string());
    if (options_.verbosity_ >= 3) {
      for (auto u : collateral) {
        logger.log(3, " with {}", u.term->to_string());
      }
    }

    // Most IC3 implementations will have only a single element in the vector
    // e.g. a single clause. But this is not guaranteed for all
    // for example, interpolant-based generalization for bit-vectors is not
    // always a single clause
    size_t min_idx = frames_.size();
    for (auto bu : collateral) {
      // try to push
      size_t idx = find_highest_frame(i, bu);
      constrain_frame(idx, bu);
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
    // collateral is a vector of predecessors
    // for now, assume there is only one
    // TODO: extend this to support multiple predecessors
    assert(collateral.size() == 1);
    add_proof_goal(collateral.at(0), i - 1, make_shared<ProofGoal>(pg));
    return false;
  }
}

bool IC3Base::is_blocked(const ProofGoal & pg)
{
  assert(solver_context_ == 0);

  push_solver_context();
  assert_frame_labels(pg.idx);
  solver_->assert_formula(pg.target.term);
  Result r = solver_->check_sat();
  pop_solver_context();

  return r.is_unsat();
}

bool IC3Base::propagate(size_t i)
{
  assert(i + 1 < frames_.size());

  unordered_set<size_t> indices_to_remove;
  const vector<IC3Formula> & Fi = frames_.at(i);

  push_solver_context();
  assert_frame_labels(i);
  assert_trans_label();

  for (size_t j = 0; j < Fi.size(); ++j) {
    const Term & t = Fi.at(j).term;

    // Relative inductiveness check
    // Check F[i] /\ t /\ T /\ -t'
    // NOTE: asserting t is redundant because t \in F[i]
    push_solver_context();
    solver_->assert_formula(solver_->make_term(Not, ts_->next(t)));

    Result r = solver_->check_sat();
    assert(!r.is_unknown());
    if (r.is_unsat()) {
      // mark for removal
      indices_to_remove.insert(j);
    }

    pop_solver_context();
  }

  pop_solver_context();

  // keep invariant that terms are kept at highest frame
  // where they are known to hold
  // need to remove some terms from the current frame
  vector<IC3Formula> new_frame_i;
  new_frame_i.reserve(Fi.size() - indices_to_remove.size());
  for (size_t j = 0; j < Fi.size(); ++j) {
    if (indices_to_remove.find(j) == indices_to_remove.end()) {
      new_frame_i.push_back(Fi.at(j));
    } else {
      // add to next frame
      constrain_frame(i + 1, Fi.at(j));
    }
  }

  frames_[i] = new_frame_i;

  return new_frame_i.empty();
}

void IC3Base::push_frame()
{
  assert(frame_labels_.size() == frames_.size());
  // pushes an empty frame
  frame_labels_.push_back(
      solver_->make_symbol("__frame_label_" + std::to_string(frames_.size()),
                           solver_->make_sort(BOOL)));
  frames_.push_back({});
}

void IC3Base::constrain_frame(size_t i, const IC3Formula & constraint)
{
  assert(solver_context_ == 0);
  assert(i > 0);  // there's a special case for frame 0
  assert(i < frame_labels_.size());
  assert(frame_labels_.size() == frames_.size());
  solver_->assert_formula(
      solver_->make_term(Implies, frame_labels_.at(i), constraint.term));
  frames_.at(i).push_back(constraint);
}

void IC3Base::assert_frame_labels(size_t i) const
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

Term IC3Base::get_frame_term(size_t i) const
{
  // TODO: decide if frames should hold IC3Formulas or terms
  //       need to special case initial state if using IC3Formulas
  if (i == 0) {
    // F[0] is always the initial states constraint
    return ts_->init();
  }

  Term res = solver_true_;
  for (size_t j = i; j < frames_.size(); ++j) {
    for (auto u : frames_[j]) {
      res = solver_->make_term(And, res, u.term);
    }
  }
  return res;
}

void IC3Base::assert_trans_label() const
{
  // shouldn't be a scenario where trans is asserted at base context
  // just because of how IC3 works
  assert(solver_context_ > 0);
  solver_->assert_formula(trans_label_);
}

bool IC3Base::has_proof_goals() const { return !proof_goals_.empty(); }

ProofGoal IC3Base::get_next_proof_goal()
{
  assert(has_proof_goals());
  ProofGoal pg = proof_goals_.back();
  proof_goals_.pop_back();
  return pg;
}

void IC3Base::add_proof_goal(const IC3Formula & c,
                             size_t i,
                             shared_ptr<ProofGoal> n)
{
  // IC3Formula aligned with frame so proof goal should be negated
  // e.g. for bit-level IC3, IC3Formula is a Clause and the proof
  // goal should be a Cube
  assert(!c.is_disjunction());
  assert(ic3formula_check_valid(c));
  proof_goals_.push_back(ProofGoal(c, i, n));
}

bool IC3Base::check_intersects(const Term & A, const Term & B)
{
  // should only do this check starting from context 0
  // don't want polluting assumptions
  assert(solver_context_ == 0);
  push_solver_context();
  solver_->assert_formula(A);
  solver_->assert_formula(B);
  Result r = solver_->check_sat();
  pop_solver_context();
  return r.is_sat();
}

bool IC3Base::check_intersects_initial(const Term & t)
{
  return check_intersects(init_label_, t);
}

void IC3Base::fix_if_intersects_initial(TermVec & to_keep, const TermVec & rem)
{
  if (rem.size() != 0) {
    // TODO: there's a tricky issue here. The reducer doesn't have the label
    // assumptions so we can't use init_label_ here. need to come up with a
    // better interface. Should we add label assumptions to reducer?
    Term formula = solver_->make_term(And, ts_->init(), make_and(to_keep));
    reducer_.reduce_assump_unsatcore(formula,
                                     rem,
                                     to_keep,
                                     NULL,
                                     options_.ic3_gen_max_iter_,
                                     options_.random_seed_);
  }
}

size_t IC3Base::find_highest_frame(size_t i, const IC3Formula & u)
{
  assert(u.is_disjunction());
  Term c = u.term;
  push_solver_context();
  solver_->assert_formula(c);
  solver_->assert_formula(solver_->make_term(Not, ts_->next(c)));
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

Term IC3Base::make_and(TermVec vec, SmtSolver slv) const
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

void IC3Base::push_solver_context()
{
  solver_->push();
  solver_context_++;
}

void IC3Base::pop_solver_context()
{
  solver_->pop();
  solver_context_--;
}

Term IC3Base::label(const Term & t)
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

smt::Term IC3Base::smart_not(const Term & t) const
{
  Op op = t->get_op();
  if (op == Not) {
    TermVec children(t->begin(), t->end());
    assert(children.size() == 1);
    return children[0];
  } else {
    return solver_->make_term(Not, t);
  }
}

}  // namespace pono
