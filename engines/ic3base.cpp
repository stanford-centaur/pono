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

/**
 * Priority queue of proof obligations inspired by open-source ic3ia
 * implementation
 */
class ProofGoalQueue
{
 public:
  ~ProofGoalQueue() { clear(); }

  void clear()
  {
    for (auto p : store_) {
      delete p;
    }
    store_.clear();
    while (!queue_.empty()) {
      queue_.pop();
    }
  }

  void new_proof_goal(const IC3Formula & c,
                      unsigned int t,
                      const ProofGoal * n = NULL)
  {
    ProofGoal * pg = new ProofGoal(c, t, n);
    queue_.push(pg);
    store_.push_back(pg);
  }

  ProofGoal * top() { return queue_.top(); }
  void pop() { queue_.pop(); }
  bool empty() const { return queue_.empty(); }

 private:
  typedef std::
      priority_queue<ProofGoal *, std::vector<ProofGoal *>, ProofGoalOrder>
          Queue;
  Queue queue_;
  std::vector<ProofGoal *> store_;
};

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

/** Syntactic subsumption check for clauses: ? a subsumes b ?
 *  @param IC3Formula a
 *  @param IC3Formula b
 *  returns true iff 'a subsumes b'
 */
static bool subsumes(const IC3Formula &a, const IC3Formula &b)
{
  assert(a.disjunction);
  assert(a.disjunction == b.disjunction);
  const TermVec &ac = a.children;
  const TermVec &bc = b.children;
  // NOTE: IC3Formula children are sorted on construction
  //       Uses unique id of term, from term->get_id()
  return ac.size() <= bc.size()
         && std::includes(bc.begin(), bc.end(), ac.begin(), ac.end());
}

/** IC3Base */

IC3Base::IC3Base(const Property & p, const TransitionSystem & ts,
                 const SmtSolver & s, PonoOptions opt)
    : super(p, ts, s, opt),
      reducer_(create_solver(s->get_solver_enum())),
      solver_context_(0),
      num_check_sat_since_reset_(0),
      failed_to_reset_solver_(false),
      cex_pg_(nullptr)
{
}

IC3Base::~IC3Base()
{
  if (cex_pg_) {
    delete cex_pg_;
    cex_pg_ = nullptr;
  }
}

void IC3Base::initialize()
{
  if (initialized_) {
    return;
  }

  // abstract the transition relation if this is a CEGAR implementation
  // otherwise it is a No-Op
  abstract();

  super::initialize();

  // check whether this flavor of IC3 can be applied to this transition system
  check_ts();

  assert(solver_context_ == 0);  // expecting to be at base context level
  solver_true_ = solver_->make_term(true);

  frames_.clear();
  frame_labels_.clear();
  // first frame is always the initial states
  push_frame();
  // can't use constrain_frame for initial states because not guaranteed to be
  // an IC3Formula it's handled specially
  solver_->assert_formula(
      solver_->make_term(Implies, frame_labels_.at(0), ts_.init()));
  push_frame();

  // set semantics of TS labels
  Sort boolsort = solver_->make_sort(BOOL);
  assert(!init_label_);
  assert(!trans_label_);
  // frame 0 label is identical to init label
  init_label_ = frame_labels_[0];
  trans_label_ = solver_->make_symbol("__trans_label", boolsort);
  solver_->assert_formula(
      solver_->make_term(Implies, trans_label_, ts_.trans()));
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
  assert(reached_k_ + 1 >= 0);
  for (size_t i = reached_k_ + 1; i <= k; ++i) {
    // reset cex_pg_ to null
    // there might be multiple abstract traces if there's a derived class
    // doing abstraction refinement
    if (cex_pg_) {
      delete cex_pg_;
      cex_pg_ = nullptr;
    }

    res = step(i);
    if (res != ProverResult::UNKNOWN) {
      return res;
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
  return IC3Formula(term, c, true);
}

IC3Formula IC3Base::ic3formula_conjunction(const TermVec & c) const
{
  assert(c.size());
  Term term = c.at(0);
  for (size_t i = 1; i < c.size(); ++i) {
    term = solver_->make_term(And, term, c[i]);
  }
  return IC3Formula(term, c, false);
}

IC3Formula IC3Base::ic3formula_negate(const IC3Formula & u) const
{
  const TermVec & children = u.children;
  assert(!u.is_null());
  assert(children.size());

  TermVec neg_children;
  neg_children.reserve(children.size());
  Term nc = smart_not(children.at(0));

  bool is_clause = u.disjunction;
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
  return IC3Formula(term, neg_children, !is_clause);
}

bool IC3Base::intersects_bad(IC3Formula & out)
{
  push_solver_context();
  // assert the last frame (conjunction over clauses)
  assert_frame_labels(reached_k_ + 1);
  // see if it intersects with bad
  solver_->assert_formula(bad_);
  Result r = check_sat();

  if (r.is_sat()) {
    const IC3Formula &c = get_model_ic3formula();
    // reduce c
    TermVec red_c;
    reducer_.reduce_assump_unsatcore(smart_not(bad_), c.children, red_c);

    out = ic3formula_conjunction(red_c);
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
  if (!block_all()) {
    // counter-example
    return ProverResult::FALSE;
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
  Result r = check_sat();
  if (r.is_sat()) {
    const IC3Formula &c = get_model_ic3formula();
    cex_pg_ = new ProofGoal(c, 0, nullptr);
    pop_solver_context();
    return ProverResult::FALSE;
  } else {
    assert(r.is_unsat());
    reached_k_ = 0;  // keep reached_k_ aligned with number of frames
  }
  pop_solver_context();
  return ProverResult::UNKNOWN;
}

bool IC3Base::rel_ind_check(size_t i, const IC3Formula & c, IC3Formula & out)
{
  assert(i > 0);
  assert(i < frames_.size());
  // expecting to be the polarity for proof goals, not frames
  // e.g. a conjunction
  assert(!c.disjunction);

  assert(solver_context_ == 0);
  push_solver_context();

  // F[i-1]
  assert_frame_labels(i - 1);
  // -c
  solver_->assert_formula(solver_->make_term(Not, c.term));
  // Trans
  assert_trans_label();
  // c'
  solver_->assert_formula(ts_.next(c.term));

  Result r = check_sat();
  if (r.is_sat()) {
    if (options_.ic3_pregen_) {
      out = generalize_predecessor(i, c);
    } else {
      out = get_model_ic3formula();
    }
    assert(ic3formula_check_valid(out));
    pop_solver_context();
  } else {
    // TODO: consider automatically taking advantage
    //       of an unsat core. Took it out for now (was in MBIC3)
    //       because it needs to work for any IC3Formula
    //       Maybe IC3Formula needs to know how to generalize itself
    //         or at least how to make a conjunctive partition
    //         or it's possible they all can function approximately the same
    //       would also have to move the pop_solver_context later
    out = c;
    pop_solver_context();
  }
  assert(!solver_context_);

  if (r.is_sat()) {
    // this check needs to be here after the solver context has been popped
    // if i == 1 and there's a predecessor, then it should be an initial state
    assert(i != 1 || check_intersects_initial(out.term));

    // should never intersect with a frame before F[i-1]
    // otherwise, this predecessor should have been found
    // in a previous step (before a new frame was pushed)
    assert(i < 2 || !check_intersects(out.term, get_frame_term(i - 2)));
  }

  assert(!r.is_unknown());
  return r.is_unsat();
}

// Helper methods

bool IC3Base::block_all()
{
  assert(!solver_context_);
  ProofGoalQueue proof_goals;
  IC3Formula bad_goal;
  while (intersects_bad(bad_goal)) {
    assert(bad_goal.term);  // expecting non-null
    assert(proof_goals.empty());  // bad should be the first goal each iteration
    proof_goals.new_proof_goal(bad_goal, frontier_idx(), nullptr);

    while (!proof_goals.empty()) {
      const ProofGoal * pg = proof_goals.top();

      if (!pg->idx) {
        // went all the way back to initial
        // TODO refactor refinement to not use cex_pg_
        // need to create a new proof goal that's not managed by the queue
        cex_pg_ = new ProofGoal(pg->target, pg->idx, pg->next);
        RefineResult s = refine();
        if (s == REFINE_SUCCESS) {
          // on successful refinement, clear the queue of proof goals
          // which might not have been precise
          // TODO might have to change this if there's an algorithm
          // that refines but can keep proof goals around
          proof_goals.clear();

          // and reset cex_pg_
          if (cex_pg_) {
            delete cex_pg_;
            cex_pg_ = nullptr;
          }
          continue;
        } else if (s == REFINE_NONE) {
          // this is a real counterexample
          // TODO refactor this
          assert(cex_pg_);
          assert(cex_pg_->target.term == pg->target.term);
          assert(cex_pg_->idx == pg->idx);
          return false;
        } else {
          assert(s == REFINE_FAIL);
          throw PonoException("Refinement failed");
        }
      }

      if (is_blocked(pg)) {
        logger.log(3,
                   "Skipping already blocked proof goal <{}, {}>",
                   pg->target.term,
                   pg->idx);
        // remove the proof goal since it has already been blocked
        assert(pg == proof_goals.top());
        proof_goals.pop();
        continue;
      }

      IC3Formula collateral;  // populated by rel_ind_check
      if (rel_ind_check(pg->idx, pg->target, collateral)) {
        // this proof goal can be blocked
        assert(!solver_context_);
        assert(collateral.term);
        logger.log(
            3, "Blocking term at frame {}: {}", pg->idx, pg->target.term);

        // remove the proof goal now that it has been blocked
        assert(pg == proof_goals.top());
        proof_goals.pop();

        if (options_.ic3_indgen_) {
          collateral = inductive_generalization(pg->idx, collateral);
        } else {
          // just negate the term
          collateral = ic3formula_negate(collateral);
        }

        size_t idx = find_highest_frame(pg->idx, collateral);
        assert(idx >= pg->idx);
        constrain_frame(idx, collateral);

        // re-add the proof goal at a higher frame if not blocked
        // up to the frontier
        if (idx < frontier_idx()) {
          assert(!pg->target.disjunction);
          proof_goals.new_proof_goal(pg->target, idx + 1, pg->next);
        }

      } else {
        // could not block this proof goal
        assert(collateral.term);
        proof_goals.new_proof_goal(collateral, pg->idx - 1, pg);
      }
    }  // end while(!proof_goals.empty())

    assert(!(bad_goal = IC3Formula()).term);  // in debug mode, reset it
  }                                           // end while(intersects_bad())

  assert(proof_goals.empty());
  return true;
}

bool IC3Base::is_blocked(const ProofGoal * pg)
{
  // syntactic check
  for (size_t i = pg->idx; i < frames_.size(); ++i) {
    const vector<IC3Formula> & Fi = frames_.at(i);
    for (size_t j = 0; j < Fi.size(); ++j) {
      if (subsumes(Fi[j], ic3formula_negate(pg->target))) {
        return true;
      }
    }
  }

  // now semantic check
  assert(solver_context_ == 0);

  push_solver_context();
  assert_frame_labels(pg->idx);
  solver_->assert_formula(pg->target.term);
  Result r = check_sat();
  pop_solver_context();

  return r.is_unsat();
}

bool IC3Base::propagate(size_t i)
{
  assert(i + 1 < frames_.size());

  vector<IC3Formula> to_push;
  vector<IC3Formula> & Fi = frames_.at(i);

  push_solver_context();
  assert_frame_labels(i);
  assert_trans_label();

  size_t k = 0;
  for (size_t j = 0; j < Fi.size(); ++j) {
    const Term & t = Fi.at(j).term;

    // Relative inductiveness check
    // Check F[i] /\ t /\ T /\ -t'
    // NOTE: asserting t is redundant because t \in F[i]
    push_solver_context();
    solver_->assert_formula(solver_->make_term(Not, ts_.next(t)));

    Result r = check_sat();
    assert(!r.is_unknown());
    if (r.is_unsat()) {
      to_push.push_back(Fi.at(j));
    } else {
      Fi[k++] = Fi.at(j);
    }

    pop_solver_context();
  }
  Fi.resize(k);

  pop_solver_context();

  for (const auto &f : to_push) {
    constrain_frame(i + 1, f, false);
  }

  return Fi.empty();
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

void IC3Base::constrain_frame(size_t i, const IC3Formula & constraint,
                              bool new_constraint)
{
  assert(solver_context_ == 0);
  assert(i < frame_labels_.size());
  assert(constraint.disjunction);
  assert(ts_.only_curr(constraint.term));

  if (new_constraint) {
    for (size_t j = 1; j <= i; ++j) {
      vector<IC3Formula> & Fj = frames_.at(j);
      size_t k = 0;
      for (size_t l = 0; l < Fj.size(); ++l) {
        if (!subsumes(constraint, Fj[l])) {
          Fj[k++] = Fj[l];
        }
      }
      Fj.resize(k);
    }
  }

  assert(i > 0);  // there's a special case for frame 0

  constrain_frame_label(i, constraint);
  frames_.at(i).push_back(constraint);
}

void IC3Base::constrain_frame_label(size_t i, const IC3Formula & constraint)
{
  assert(frame_labels_.size() == frames_.size());

  solver_->assert_formula(
      solver_->make_term(Implies, frame_labels_.at(i), constraint.term));
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
    return ts_.init();
  }

  Term res = solver_true_;
  for (size_t j = i; j < frames_.size(); ++j) {
    for (const auto &u : frames_[j]) {
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

bool IC3Base::check_intersects(const Term & A, const Term & B)
{
  // should only do this check starting from context 0
  // don't want polluting assumptions
  assert(solver_context_ == 0);
  push_solver_context();
  solver_->assert_formula(A);
  solver_->assert_formula(B);
  Result r = check_sat();
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
    const Term &formula = solver_->make_term(And, ts_.init(), make_and(to_keep));
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
  assert(u.disjunction);
  const Term &c = u.term;
  push_solver_context();
  solver_->assert_formula(c);
  solver_->assert_formula(solver_->make_term(Not, ts_.next(c)));
  assert_trans_label();

  Result r;
  size_t j = i;
  for (; j + 1 < frames_.size(); ++j) {
    push_solver_context();
    assert_frame_labels(j);
    r = check_sat();
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

void IC3Base::reset_solver()
{
  assert(solver_context_ == 0);

  if (failed_to_reset_solver_) {
    // don't even bother trying
    // this solver doesn't support reset_assertions
    return;
  }

  try {
    solver_->reset_assertions();

    // Now need to add back in constraints at context level 0
    logger.log(2, "IC3Base: Reset solver and now re-adding constraints.");

    // define init and trans label
    assert(init_label_ == frame_labels_.at(0));
    solver_->assert_formula(
        solver_->make_term(Implies, init_label_, ts_.init()));
    solver_->assert_formula(
        solver_->make_term(Implies, trans_label_, ts_.trans()));

    for (size_t i = 0; i < frames_.size(); ++i) {
      for (const auto & constraint : frames_.at(i)) {
        constrain_frame_label(i, constraint);
      }
    }
  }
  catch (SmtException & e) {
    logger.log(1,
               "Failed to reset solver (underlying solver must not support "
               "it). Disabling solver resets for rest of run.");
    failed_to_reset_solver_ = true;
  }

  num_check_sat_since_reset_ = 0;
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
  const Op &op = t->get_op();
  if (op == Not) {
    TermVec children(t->begin(), t->end());
    assert(children.size() == 1);
    return children[0];
  } else {
    return solver_->make_term(Not, t);
  }
}

}  // namespace pono
