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
#include <numeric>

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

IC3Base::IC3Base(Property & p, const SmtSolver & s, PonoOptions opt)
    : super(p, s, opt),
      reducer_(create_solver(s->get_solver_enum())),
      solver_context_(0),
      num_check_sat_since_reset_(0),
      failed_to_reset_solver_(false),
      cex_pg_(nullptr),
      rng_(options_.random_seed_)
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
  frames_.push_back({});
  frame_labels_.push_back(
      solver_->make_symbol("__frame_label_0", solver_->make_sort(BOOL)));
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
  bad_label_ = solver_->make_symbol("__bad_label", boolsort);
  solver_->assert_formula(solver_->make_term(Implies, bad_label_, bad_));
}

ProverResult IC3Base::check_until(int k)
{
  initialize();
  // make sure derived class implemented initialize and called
  // this version of initialize with super::initialize or
  // (for experts only) set the initialized_ flag without
  // ever initializing base classes
  assert(initialized_);

  if (step_0() == ProverResult::FALSE) {
    return ProverResult::FALSE;
  }

  while (depth() <= k) {
    IC3Formula bad;
    while (get_bad(bad)) {
      RefineResult s = rec_block(bad);
      if (s == REFINE_NONE) {
        logger.log(1, "found counterexample at depth {}", depth());
        return ProverResult::FALSE;
      } else if (s == REFINE_FAIL) {
        logger.log(1, "unknown result at depth {}", depth());
        return ProverResult::UNKNOWN;
      }
    }
    reached_k_ = depth();
    push_frame();
    if (propagate()) {
      return ProverResult::TRUE;
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
  return IC3Formula(term, neg_children, !is_clause);
}

bool IC3Base::get_bad(IC3Formula & out)
{
  push_solver_context();

  // assert the last frame (conjunction over clauses)
  assert_frame_labels(reached_k_ + 1);
  // see if it intersects with bad
  activate_bad();

  Result r = check_sat();
  if (r.is_sat()) {
    out = get_model_ic3formula();
    generalize_bad(out);
  }

  pop_solver_context();

  assert(!r.is_unknown());
  return r.is_sat();
}

ProverResult IC3Base::step(int i)
{
  throw PonoException("IC3Base::step removed");
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

// Helper methods

RefineResult IC3Base::rec_block(const IC3Formula & bad)
{
  ProofGoalQueue queue;
  queue.push_new(bad, depth());

  while (!queue.empty()) {
    // periodically reset the solver -- clean up "garbage"
    // (e.g. subsumed/learned clauses, bad variable scores...)
    if (options_.ic3_reset_interval_
        && num_check_sat_since_reset_ >= options_.ic3_reset_interval_) {
      reset_solver();
    }

    ProofGoal * p = queue.top();

    logger.log(3,
               "looking at proof goal of size {} at idx {}",
               p->target.children.size(),
               p->idx);
    logger.log(4, "{}", p->target.term);

    if (p->idx == 0) {
      // doing differently than msat-ic3ia -- reconstructing in refine()
      cex_pg_ = p;
      RefineResult s = refine();
      if (s == REFINE_SUCCESS) {
        // upon successful refinement, we clear the queue of proof
        // obligations. This is because we have added more predicates,
        // so the proof obligations still in the queue now might be
        // imprecise wrt. the current predicate abstraction. If we
        // keep them around, we might get spurious counterexamples
        // even if the predicate abstraction is precise enough. In
        // principle we could handle this, but it is simpler to just
        // flush the queue
        while (!queue.empty()) {
          queue.pop();
        }
      }
      return s;
    }

    assert(!check_intersects_initial(p->target.term));

    if (!is_blocked(p)) {
      IC3Formula c;
      if (block(p->target, p->idx, &c, true)) {
        logger.log(3,
                   "CTI successfully blocked: subcube of size {} is inductive "
                   "relative to frame {}",
                   c.children.size(),
                   p->idx - 1);
        logger.log(4, "{}", c.term);

        // c is inductive relative to F[idx-1]
        unsigned int idx = p->idx;
        // generalize c and find the highest idx frame for which
        // relative induction holds
        generalize_and_push(c, idx);
        // add c to F[idx]...F[0]
        constrain_frame(c, idx);
        if (idx < depth()) {
          // if we are not at the frontier, try to block p at later
          // steps. Remember that p is a state that leads to a bad
          // state in depth()-p->idx steps. Therefore, p is also
          // bad and if we want to prove the property, our inductive
          // invariant must not intersect with it. Since we have p
          // already, it makes sense to exploit this information and
          // schedule the proof obligation right away. Notice that
          // this just an optimization, and it can be turned off
          // without affecting correctness or completeness. However,
          // this optimization is typically quite effective, and it
          // is the reason why IC3 might sometimes find cex traces
          // that are longer that the current depth()
          queue.push_new(p->target, p->idx + 1, p->next);
        }
        queue.pop();
      } else {
        logger.log(3, "got CTI of size {}:", c.children.size());
        logger.log(4, "{}", c.term);
        // c is a predecessor to the bad cube p->cube, so we need to
        // block it at earlier steps in the trace
        queue.push_new(c, p->idx - 1, p);
      }
    } else {
      queue.pop();
    }
  }

  return REFINE_SUCCESS;
}

bool IC3Base::block(const IC3Formula & c,
                    unsigned int idx,
                    IC3Formula * out,
                    bool compute_cti)
{
  assert(idx > 0);

  // check whether ~c is inductive relative to F[idx-1], i.e.
  // ~c & F[idx-1] & T |= ~c', that is
  // solve(~c & F[idx-1] & T & c') is unsat

  // activate T and F[idx-1]
  assert_frame_labels(idx - 1);
  assert_trans_label();

  assert(solver_context_ == 0);
  // assume c'
  TermVec assumps;
  IC3Formula primed = get_next(c);
  if (options_.random_seed_) {
    std::vector<size_t> idx(primed.children.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::shuffle(idx.begin(), idx.end(), rng_);

    for (size_t i : idx) {
      assumps.push_back(label(primed.children.at(i)));
    }
  } else {
    for (const auto & l : primed.children) {
      assumps.push_back(label(l));
    }
  }

  // temporarily assert ~c
  push_solver_context();
  solver_->assert_formula(ic3formula_negate(c).term);
  Result r = solver_->check_sat_assuming(assumps);
  if (r.is_unsat()) {
    // relative induction succeeds. If required (out != NULL), generalize
    // ~c to a stronger clause, by looking at the literals of c' that occur
    // in the unsat core. If g' is a subset of c',
    // then if "~c & F[idx-1] & T & g'" is unsat, then so is
    // "~g & F[idx-1] & T & g'" (since ~g is stronger than ~c)
    if (out) {
      // try minimizing using the unsat core
      UnorderedTermSet core;
      solver_->get_unsat_core(core);
      TermVec & candidate = out->children;
      TermVec rest;
      candidate.clear();
      for (size_t i = 0; i < primed.children.size(); ++i) {
        Term a = label(primed.children[i]);
        if (core.find(a) != core.end()) {
          candidate.push_back(c.children.at(i));
        } else {
          rest.push_back(c.children.at(i));
        }
      }
      // now candidate is a subset of c that is enough for
      // "~c & F[idx-1] & T & candidate'" to be unsat.
      // However, we are not done yet. If candidate intersects the
      // initial states (i.e. "init & candidate" is sat), then ~candidate
      // is not inductive relative to F[idx-1], as the base case
      // fails. We fix this by re-adding back literals until
      // "init & candidate" is unsat
      fix_if_intersects_initial(candidate, rest);
      *out = ic3formula_conjunction(candidate);
      pop_solver_context();
    } else {
      pop_solver_context();
    }
    return true;
  } else {
    // relative induction fails. If requested, extract a predecessor of c
    // (i.e. a counterexample to induction - CTI) from the model found by
    // the SMT solver
    TermVec inputs;
    if (compute_cti) {
      assert(out);
      *out = get_model_ic3formula(&inputs);
    }
    pop_solver_context();
    // TODO check that msat-ic3ia really doesn't generalize
    // doesn't generalize by default
    // if (compute_cti) {
    //   generalize_pre(primed, inputs, *out);
    // }

    assert(solver_context_ == 0);
    return false;
  }
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

bool IC3Base::propagate()
{
  std::vector<IC3Formula> to_add;

  size_t k = 1;
  for (; k < depth(); ++k) {
    to_add.clear();
    vector<IC3Formula> & f = frames_[k];
    // forward propagation: try to see if f[i] is inductive relative to
    // F[k+1]
    for (size_t i = 0; i < f.size(); ++i) {
      to_add.push_back(IC3Formula());

      logger.log(
          3, "trying to propagate cube of size {}", f.at(i).children.size());
      logger.log(4, "{}", f.at(i).term);
      logger.log(3, "from {} to {}", k, k + 1);

      if (!block(f[i], k + 1, &to_add.back(), false)) {
        to_add.pop_back();
      } else {
        logger.log(3, "success");
      }
    }

    for (IC3Formula & c : to_add) {
      constrain_frame(c, k + 1);
    }
    if (frames_[k].empty()) {
      // fixpoint: frames_[k] == frames_[k+1]
      break;
    }
  }

  if (k < depth()) {
    logger.log(2, "fixpoint found at frame {}", k);
    assert(!frames_.at(k).size());
    invar_ = get_frame_term(k + 1);
    logger.log(2, "invariant: {}", invar_);
    return true;
  }

  return false;
}

void IC3Base::push_frame()
{
  size_t depth = frames_.size() - 1;

  if (depth) {
    string row;
    for (size_t i = 1; i <= depth; ++i) {
      row += " " + std::to_string(frames_[i].size());
    }
    logger.log(1, "{} : {}", depth, row);
  }

  assert(frame_labels_.size() == frames_.size());
  // pushes an empty frame
  frame_labels_.push_back(
      solver_->make_symbol("__frame_label_" + std::to_string(frames_.size()),
                           solver_->make_sort(BOOL)));
  frames_.push_back({});
}

void IC3Base::constrain_frame(const IC3Formula & c, size_t idx)
{
  assert(solver_context_ == 0);
  assert(idx < frame_labels_.size());

  // copied from msat-ic3ia (as several other functions in this branch were)
  // trying to debug performance

  // whenever we add a clause ~c to an element of F, we also remove subsumed
  // clauses. This automatically keeps frames_ in a "delta encoded" form, in
  // which each clause is stored only in the last frame in which it
  // occurs. However, this does not remove subsumed clauses from the
  // underlying SMT solver. We address this by resetting the solver every
  // once in a while (see the comment in rec_block())
  for (size_t d = 1; d < idx + 1; ++d) {
    vector<IC3Formula> & fd = frames_[d];
    size_t j = 0;
    for (size_t i = 0; i < fd.size(); ++i) {
      if (!subsumes(c, fd[i])) {
        fd[j++] = fd[i];
      }
    }
    fd.resize(j);
  }

  solver_->assert_formula(
      solver_->make_term(Implies, frame_labels_.at(idx), c.term));
  frames_[idx].push_back(c);

  logger.log(3, "adding cube of size {} at level {}", c.children.size(), idx);
  logger.log(4, "{}", c.term);
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
    return ts_->init();
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

void IC3Base::add_proof_goal(const IC3Formula & c,
                             size_t i,
                             const ProofGoal * n)
{
  // IC3Formula aligned with frame so proof goal should be negated
  // e.g. for bit-level IC3, IC3Formula is a Clause and the proof
  // goal should be a Cube
  assert(!c.is_disjunction());
  assert(ic3formula_check_valid(c));
  proof_goals_.push_new(c, i, n);
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
    const Term &formula = solver_->make_term(And, ts_->init(), make_and(to_keep));
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
  const Term &c = u.term;
  push_solver_context();
  solver_->assert_formula(c);
  solver_->assert_formula(solver_->make_term(Not, ts_->next(c)));
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
        solver_->make_term(Implies, init_label_, ts_->init()));
    solver_->assert_formula(
        solver_->make_term(Implies, trans_label_, ts_->trans()));

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

// temporarily added from msat-ic3ia

void IC3Base::generalize(IC3Formula & c, unsigned int & idx)
{
  tmp_ = c.children;
  IC3Formula tmp_form;
  gen_needed_.clear();

  typedef std::uniform_int_distribution<int> RandInt;
  RandInt dis;

  logger.log(3,
             "trying to generalize cube of size {} at {}: ",
             c.children.size(),
             idx);
  logger.log(4, "{}", c.term);

  // ~c is inductive relative to idx-1, and we want to generalize ~c to a
  // stronger clause ~g. We do this by dropping literals and calling block()
  // again: every time block() succeeds, it will further generalize its
  // input cube using the unsat core found by the SMT solver (see above)
  //
  // More sophisticated (more effective but harder to understand) strategies
  // exist, see e.g. the paper
  // - Hassan, Somenzi, Bradley: Better generalization in IC3. FMCAD'13
  //
  for (size_t i = 0; i < tmp_.size() && tmp_.size() > 1;) {
    // randomly pick the next literal to drop
    size_t j = (options_.random_seed_
                    ? dis(rng_, RandInt::param_type(1, tmp_.size())) - 1
                    : i);
    Term l = tmp_.at(j);
    if (gen_needed_.find(l) == gen_needed_.end()) {
      auto it = tmp_.erase(tmp_.begin() + j);

      logger.log(4, "trying to drop {}", l);

      // after editing tmp, recreate it
      tmp_form = ic3formula_conjunction(tmp_);
      if (check_intersects_initial(tmp_form.term)
          || !block(tmp_form, idx, &tmp_form, false)) {
        // remember that we failed to remove l, so that we do not try
        // this again later in the loop
        gen_needed_.insert(l);
        tmp_ = tmp_form.children;
        tmp_.insert(it, l);
        ++i;
      }
    }
  }

  // after editing tmp, recreate it
  tmp_form = ic3formula_conjunction(tmp_);
  std::swap(tmp_form, c);
}

void IC3Base::push(IC3Formula & c, unsigned int & idx)
{
  IC3Formula tmp;
  // find the highest idx frame in F which can successfully block c. As a
  // byproduct, this also further strenghens ~c if possible
  while (idx < depth() - 1) {
    tmp = IC3Formula();
    if (block(c, idx + 1, &tmp, false)) {
      std::swap(tmp, c);
      ++idx;
    } else {
      break;
    }
  }
}

void IC3Base::generalize_and_push(IC3Formula & c, unsigned int & idx)
{
  generalize(c, idx);
  push(c, idx);
}

inline void IC3Base::generalize_bad(IC3Formula & c)
{
  UnorderedTermMap label_to_term;
  TermVec assumps;
  for (const auto & l : c.children) {
    Term lbl = label(l);
    label_to_term[lbl] = l;
    assumps.push_back(lbl);
  }
  push_solver_context();
  Term prop = property_.prop();
  solver_->assert_formula(prop);
  Result r = solver_->check_sat_assuming(assumps);
  if (r.is_unsat()) {
    UnorderedTermSet core;
    solver_->get_unsat_core(core);
    auto it =
        std::remove_if(assumps.begin(), assumps.end(), [&core](const Term & l) {
          return core.find(l) == core.end();
        });
    assumps.resize(it - assumps.begin());
    TermVec cube_lits;
    for (const auto & lbl : assumps) {
      cube_lits.push_back(label_to_term.at(lbl));
    }
    c = ic3formula_conjunction(cube_lits);
  }
  pop_solver_context();
}

}  // namespace pono
