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
#include "utils/term_analysis.h"

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

IC3Base::IC3Base(const Property & p,
                 const TransitionSystem & ts,
                 const SmtSolver & s,
                 PonoOptions opt)
    : super(p, ts, s, opt),
      // NOTE: this is a hack
      // TODO fix this
      reducer_(create_reducer_for(
          s->get_solver_enum(), Engine::IC3IA_ENGINE, false)),
      solver_context_(0),
      num_check_sat_since_reset_(0),
      failed_to_reset_solver_(false),
      cex_pg_(nullptr),
      boolsort_(solver_->make_sort(BOOL))
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
  reducer_.assume_label(frame_labels_.at(0), ts_.init());
  push_frame();

  // set semantics of TS labels
  assert(!init_label_);
  assert(!trans_label_);
  assert(!bad_label_);
  // frame 0 label is identical to init label
  init_label_ = frame_labels_[0];

  trans_label_ = solver_->make_symbol("__trans_label", boolsort_);
  solver_->assert_formula(
      solver_->make_term(Implies, trans_label_, ts_.trans()));
  reducer_.assume_label(trans_label_, ts_.trans());

  bad_label_ = solver_->make_symbol("__bad_label", boolsort_);
  solver_->assert_formula(solver_->make_term(Implies, bad_label_, bad_));
  reducer_.assume_label(bad_label_, bad_);
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

IC3Formula IC3Base::inductive_generalization(size_t i, const IC3Formula & c)
{
  assert(!solver_context_);
  assert(i <= frontier_idx());
  assert(!c.disjunction);  // expecting a cube
  // be default will try to find a minimal cube
  // NOTE: not necessarily minimum (e.g. it's a local minimum)

  logger.log(
      3, "trying to generalize an IC3Formula of size {}", c.children.size());

  // TODO use unsat core reducer
  // TODO use ic3_gen_max_iter_ option or remove it
  //      maybe default zero could mean unbounded
  //      seems like a good compromise

  UnorderedTermSet necessary;  // populated with children we
                               // can't drop

  IC3Formula gen = c;
  IC3Formula out;
  Term dropped;
  size_t j = 0;
  while (j < gen.children.size() && gen.children.size() > 1) {
    // TODO use random_seed_ if set for shuffling
    //      order of drop attempts

    // try dropping j
    dropped = gen.children.at(j);
    if (necessary.find(dropped) != necessary.end()) {
      // can't drop this one
      j++;
      continue;
    }

    gen.children.erase(gen.children.begin() + j);

    // TODO: decide if it's too expensive to create fresh
    //       IC3Formula each time -- which sorts the elements
    //       if so, could consider not automatically sorting
    //       and instead only doing it for subsumption checks
    gen = ic3formula_conjunction(gen.children);

    if (!check_intersects_initial(gen.term)
        && rel_ind_check(i, gen, out, false)) {
      // we can drop this literal

      // out was generalized with an unsat core in
      // rel_ind_check
      // we can't rely on the order of the children
      // being the same
      gen = out;
      j = 0;  // start iteration over
    } else {
      // could not drop this child
      necessary.insert(dropped);
      // NOTE gen.term won't be updated
      //      but gen will be reconstructed in
      //      next iteration anyway
      gen.children.push_back(dropped);

      // NOTE: don't need to increment j because
      //       the one at position j was put at
      //       end of vector
      assert(j + 1 == gen.children.size() || gen.children.at(j) != dropped);
    }
  }

  // reconstruct the IC3Formula -- need to make sure term is valid
  // since we've been modifying gen.children
  gen = ic3formula_conjunction(gen.children);
  assert(!check_intersects_initial(gen.term));
  IC3Formula block = ic3formula_negate(gen);
  assert(block.disjunction);
  return block;
}

void IC3Base::predecessor_generalization(size_t i,
                                         const IC3Formula & c,
                                         IC3Formula & pred)
{
  // by default does no generalization
  return;
}

bool IC3Base::intersects_bad(IC3Formula & out)
{
  push_solver_context();
  // assert the last frame (conjunction over clauses)
  assert_frame_labels(reached_k_ + 1);
  // see if it intersects with bad
  solver_->assert_formula(bad_label_);
  // don't need transition relation for this check
  // can deactivate it
  solver_->assert_formula(solver_->make_term(Not, trans_label_));
  Result r = check_sat();

  if (r.is_sat()) {
    out = get_model_ic3formula();
    assert(out.term);
    assert(out.children.size());
    assert(ic3formula_check_valid(out));

    // reduce
    TermVec red_c;
    // with abstraction can't guarantee this is unsat
    if (reducer_.reduce_assump_unsatcore(
            smart_not(bad_), out.children, red_c)) {
      logger.log(2,
                 "generalized bad cube to {}/{}",
                 red_c.size(),
                 out.children.size());
      out = ic3formula_conjunction(red_c);

      assert(out.term);
      assert(out.children.size());
      assert(ic3formula_check_valid(out));
    } else {
      logger.log(2, "generalizing bad failed");
    }
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
  for (size_t j = 1; j < frontier_idx(); ++j) {
    if (propagate(j)) {
      assert(j + 1 < frames_.size());
      // save the invariant
      // which is the frame that just had all terms
      // from the previous frames propagated
      set_invar(j + 1);
      return ProverResult::TRUE;
    }
  }

  reset_solver();

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

bool IC3Base::rel_ind_check(size_t i,
                            const IC3Formula & c,
                            IC3Formula & out,
                            bool get_pred)
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

  // use assumptions for c' so we can get cheap initial
  // generalization if the check is unsat

  // NOTE: relying on same order between assumps_ and c.children
  assumps_.clear();
  {
    // TODO shuffle assumps and (a copy of) c.children
    //      if random seed is set
    Term lbl, ccnext;
    for (const auto & cc : c.children) {
      ccnext = ts_.next(cc);
      lbl = label(ccnext);
      if (lbl != ccnext) {
        // only need to add assertion if the label is not the same as ccnext
        // could be the same if ccnext is already a literal
        solver_->assert_formula(solver_->make_term(Implies, lbl, ccnext));
      }
      assumps_.push_back(lbl);
    }
  }

  Result r = check_sat_assuming(assumps_);
  if (r.is_sat()) {
    if (get_pred) {
      out = get_model_ic3formula();
      if (options_.ic3_pregen_) {
        predecessor_generalization(i, c, out);
        assert(out.term);
        assert(out.children.size());
        assert(!out.disjunction);  // expecting a conjunction
      }
    }
    assert(ic3formula_check_valid(out));
    pop_solver_context();
  } else {
    // Use unsat core to get cheap generalization

    UnorderedTermSet core;
    solver_->get_unsat_core(core);
    assert(core.size());

    TermVec gen;  // cheap unsat-core generalization of c
    TermVec rem;  // conjuncts removed by unsat core
                  // might need to be re-added if it
                  // ends up intersecting with initial
    assert(assumps_.size() == c.children.size());
    for (size_t i = 0; i < assumps_.size(); ++i) {
      if (core.find(assumps_.at(i)) == core.end()) {
        rem.push_back(c.children.at(i));
      } else {
        gen.push_back(c.children.at(i));
      }
    }

    pop_solver_context();

    fix_if_intersects_initial(gen, rem);
    assert(gen.size() >= core.size());

    // keep it as a conjunction for now
    out = ic3formula_conjunction(gen);
  }
  assert(!solver_context_);

  if (r.is_sat() && get_pred) {
    assert(out.term);
    assert(out.children.size());

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

        assert(collateral.disjunction);
        assert(collateral.term);
        assert(collateral.children.size());
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
  assert(!solver_context_);
  assert(i < frontier_idx());

  vector<IC3Formula> & Fi = frames_.at(i);

  size_t k = 0;
  IC3Formula gen;
  for (size_t j = 0; j < Fi.size(); ++j) {
    const IC3Formula & c = Fi.at(j);
    assert(c.disjunction);
    assert(c.term);
    assert(c.children.size());

    // NOTE: rel_ind_check works on conjunctions
    //       need to negate
    if (rel_ind_check(i + 1, ic3formula_negate(c), gen, false)) {
      // can push to next frame
      // got unsat-core based generalization
      assert(gen.term);
      assert(gen.children.size());
      constrain_frame(i + 1, ic3formula_negate(gen), false);
    } else {
      // have to keep this one at this frame
      Fi[k++] = c;
    }
  }

  // get rid of garbage at end of frame
  Fi.resize(k);

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
  reducer_.assume_label(frame_labels_.at(i), constraint.term);
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
  assert(!solver_context_);
  if (rem.size() != 0) {
    Term formula = solver_->make_term(And, init_label_, make_and(to_keep));

    bool success = reducer_.reduce_assump_unsatcore(formula,
                                                    rem,
                                                    to_keep,
                                                    NULL,
                                                    options_.ic3_gen_max_iter_,
                                                    options_.random_seed_);
    assert(success);
  }
}

size_t IC3Base::find_highest_frame(size_t i, IC3Formula & u)
{
  assert(!solver_context_);
  assert(u.disjunction);
  assert(u.term);
  assert(u.children.size());

  IC3Formula conj = ic3formula_negate(u);
  IC3Formula gen;
  size_t j = i;
  for (; j < frontier_idx(); ++j) {
    assert(!conj.disjunction);
    if (rel_ind_check(j + 1, conj, gen, false)) {
      std::swap(conj, gen);
    } else {
      break;
    }
  }
  assert(!conj.disjunction);
  assert(conj.term);
  assert(conj.children.size());

  u = ic3formula_negate(conj);
  assert(u.disjunction);
  assert(u.term);
  assert(u.children.size());
  return j;
}

TermVec IC3Base::get_input_values() const
{
  TermVec out_inputs;
  out_inputs.reserve(ts_.inputvars().size());
  for (const auto & iv : ts_.inputvars()) {
    out_inputs.push_back(solver_->make_term(Equal, iv, solver_->get_value(iv)));
  }
  return out_inputs;
}

TermVec IC3Base::get_next_state_values() const
{
  TermVec out_nexts;
  out_nexts.reserve(ts_.statevars().size());
  Term nv;
  for (const auto & sv : ts_.statevars()) {
    nv = ts_.next(sv);
    out_nexts.push_back(solver_->make_term(Equal, nv, solver_->get_value(nv)));
  }
  return out_nexts;
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
    reducer_.reset_assertions();

    // Now need to add back in constraints at context level 0
    logger.log(2, "IC3Base: Reset solver and now re-adding constraints.");

    // define init, trans, and bad labels
    assert(init_label_ == frame_labels_.at(0));
    solver_->assert_formula(
        solver_->make_term(Implies, init_label_, ts_.init()));
    reducer_.assume_label(init_label_, ts_.init());

    solver_->assert_formula(
        solver_->make_term(Implies, trans_label_, ts_.trans()));
    reducer_.assume_label(trans_label_, ts_.trans());

    solver_->assert_formula(solver_->make_term(Implies, bad_label_, bad_));
    reducer_.assume_label(bad_label_, bad_);

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

  Term l;
  if (is_lit(t, boolsort_)) {
    // this can be the label itself
    l = t;
  } else {
    unsigned i = 0;
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
  }
  assert(l);

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
