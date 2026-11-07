/*********************                                                  */
/*! \file mod_ts_prop.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Functions for modifying ts and property
**
**
**/

#include "modifiers/mod_ts_prop.h"

#include <cassert>

#include "core/rts.h"
#include "smt-switch/identity_walker.h"
#include "smt-switch/utils.h"
#include "utils/exceptions.h"
#include "utils/logger.h"
#include "utils/term_analysis.h"
#include "utils/ts_manipulation.h"

using namespace smt;
using namespace std;

namespace pono {

namespace {

// Replaces floor, and integer div and mod by nonzero integer literals, with
// fresh integer variables of a transition system, and collects the bounds
// that tie each variable to the term it replaces
class FloorAndDivEliminator : public IdentityWalker
{
 public:
  FloorAndDivEliminator(TransitionSystem & ts, bool allow_new_inputvars)
      : IdentityWalker(ts.solver(), false),
        ts_(ts),
        allow_new_inputvars_(allow_new_inputvars)
  {
  }

  const TermVec & bounds() const { return bounds_; }

 protected:
  WalkerStepResult visit_term(Term & term) override
  {
    if (preorder_ || in_cache(term)) {
      return Walker_Continue;
    }

    TermVec children;
    for (const auto & c : term) {
      Term cc;
      [[maybe_unused]] bool ok = query_cache(c, cc);
      assert(ok);
      children.push_back(cc);
    }

    Op op = term->get_op();
    Term res;
    if (op.is_null()) {
      res = term;
    } else if (op.prim_op == To_Int
               && children[0]->get_sort()->get_sort_kind() == REAL) {
      res = get_replacement_var_for(To_Int, children);
    } else if ((op.prim_op == IntDiv || op.prim_op == Mod)
               && children.size() == 2 && is_nonzero_literal(children[1])) {
      // (mod t c) is t - c * (div t c), so both share the variable for div
      Term q = get_replacement_var_for(IntDiv, children);
      res = (op.prim_op == IntDiv)
                ? q
                : solver_->make_term(Minus,
                                     children[0],
                                     solver_->make_term(Mult, children[1], q));
    } else {
      res = solver_->make_term(op, children);
    }

    save_in_cache(term, res);
    return Walker_Continue;
  }

  // whether t is an integer literal other than 0, such as 3 or (- 3). A term
  // without symbols, such as (- 3 3), can still be 0, and solvers such as cvc5
  // and Z3 do not reduce it to a literal.
  // TODO: consider extending this to other terms without symbols that are
  // nonzero, such as (+ 1 2) or (* 2 (- 3))
  bool is_nonzero_literal(const Term & t) const
  {
    Term lit = (t->get_op().prim_op == Negate) ? *t->begin() : t;
    return lit->is_value() && lit != solver_->make_term(0, lit->get_sort());
  }

  // Returns the integer variable q that stands for floor(t) if op is To_Int
  // and children is { t }, and for (div t c) if op is IntDiv and children is
  // { t, c }. The first call for a term makes q and records its bounds. Later
  // calls for the same term return the same q.
  Term get_replacement_var_for(PrimOp op, TermVec children)
  {
    assert((op == To_Int && children.size() == 1)
           || (op == IntDiv && children.size() == 2));

    // a t over next-state variables only, such as x' + 1, gets q', the
    // next-state copy of the variable q for x + 1: eliminate_floor_and_div adds
    // q's bounds for both the current and the next state. The divisor c is a
    // literal and stays as it is. A t over both, such as x + x', gets a
    // variable of its own.
    const Term & t = children[0];
    Term curr_t = ts_.curr(t);
    if (curr_t != t && ts_.only_curr(curr_t) && ts_.next(curr_t) == t) {
      children[0] = curr_t;
      return ts_.next(get_replacement_var_for(op, children));
    }

    Term key = solver_->make_term(op, children);
    auto it = vars_.find(key);
    if (it != vars_.end()) {
      return it->second;
    }

    // bounds over current-state variables only constrain every state, so q
    // must be a state variable, whose next-state copy stands for t over
    // next-state variables. Other bounds constrain only the transition that t
    // occurs in, so q can be an input variable.
    const char * role = (op == IntDiv) ? "quotient" : "floor_of_real";
    Sort intsort = solver_->make_sort(INT);
    Term q = (allow_new_inputvars_ && !ts_.only_curr(t))
                 ? ts_.make_generated_inputvar(role, intsort)
                 : ts_.make_generated_statevar(role, intsort);
    Term lower, upper;
    if (op == IntDiv) {
      const Term & c = children[1];
      Term cq = solver_->make_term(Mult, c, q);
      lower = solver_->make_term(Le, cq, t);
      upper = solver_->make_term(
          Lt, t, solver_->make_term(Plus, cq, solver_->make_term(Abs, c)));
    } else {
      Term real_q = solver_->make_term(To_Real, q);
      lower = solver_->make_term(Le, real_q, t);
      upper = solver_->make_term(
          Lt,
          t,
          solver_->make_term(
              Plus, real_q, solver_->make_term(1, t->get_sort())));
    }
    bounds_.push_back(solver_->make_term(And, lower, upper));
    vars_[key] = q;
    return q;
  }

  TransitionSystem & ts_;
  bool allow_new_inputvars_;
  UnorderedTermMap vars_;  // floor or div term to its variable
  TermVec bounds_;         // lower and upper bounds of each variable
};

}  // namespace

void eliminate_floor_and_div(TransitionSystem & ts,
                             Term & prop,
                             bool allow_new_inputvars)
{
  if (ts.is_functional()) {
    throw PonoException(
        "Eliminating div and mod requires a relational transition system");
  }
  RelationalTransitionSystem & rts =
      static_cast<RelationalTransitionSystem &>(ts);

  FloorAndDivEliminator eliminator(ts, allow_new_inputvars);
  Term init = ts.init();
  Term trans = ts.trans();
  init = eliminator.visit(init);
  trans = eliminator.visit(trans);
  prop = eliminator.visit(prop);

  const TermVec & bounds = eliminator.bounds();
  if (bounds.empty()) {
    return;
  }
  logger.log(1,
             "Replaced {} floors and integer divisions with variables",
             bounds.size());

  rts.set_init(init);
  rts.set_trans(trans);
  for (const auto & b : bounds) {
    if (ts.no_next(b)) {
      ts.add_constraint(b);
    } else {
      rts.constrain_trans(b);
    }
  }
}

TransitionSystem pseudo_init_and_prop(TransitionSystem & ts, Term & prop)
{
  logger.log(1, "Modifying init and prop");

  RelationalTransitionSystem rts =
      ts.is_functional() ? RelationalTransitionSystem(ts.solver())
                         : static_cast<RelationalTransitionSystem &>(ts);

  // make it relational if it wasn't already
  if (ts.is_functional()) {
    for (const auto & sv : ts.statevars()) {
      rts.add_statevar(sv, ts.next(sv));
    }
    for (const auto & iv : ts.inputvars()) {
      rts.add_inputvar(iv);
    }
    rts.set_init(ts.init());
    rts.set_trans(ts.trans());
  }

  assert(!rts.is_functional());

  Sort boolsort = rts.make_sort(BOOL);

  // Save original initial state constraints
  TermVec init_conjuncts;
  conjunctive_partition(rts.init(), init_conjuncts, true);

  // create a pseudo initial state
  Term pseudo_init = rts.make_generated_statevar("pseudo_init", boolsort);
  Term not_pseudo_init = rts.make_term(Not, pseudo_init);

  // guard property with it
  prop = rts.make_term(Implies, not_pseudo_init, prop);

  // create a fake transition from pseudo initial state
  // and enforce that the second state is constrained by the original
  // initial state constraints
  // NOTE: need to guard transition relation as well for correctness
  TermVec trans_conjuncts;
  conjunctive_partition(rts.trans(), trans_conjuncts, true);
  // reset transition relation
  rts.set_trans(rts.make_term(true));
  // add fake transition from pseudo_init to actual initial states
  for (const auto & ic : init_conjuncts) {
    assert(rts.only_curr(ic));
    rts.constrain_trans(rts.make_term(Implies, pseudo_init, rts.next(ic)));
  }
  // add regular transitions back in
  for (const auto & tc : trans_conjuncts) {
    rts.constrain_trans(rts.make_term(Implies, not_pseudo_init, tc));
  }

  rts.set_init(pseudo_init);
  rts.assign_next(pseudo_init, rts.make_term(false));

  // now create a property monitor
  Term new_prop = rts.make_generated_statevar("prop_monitor", boolsort);
  if (rts.only_curr(prop)) {
    rts.add_invar(rts.make_term(Equal, new_prop, prop));
  } else {
    // property monitor
    rts.assign_next(new_prop, prop);
  }
  // need to assume in the pseudo initial state
  // or else there's a trivial counterexample
  rts.constrain_init(new_prop);

  // update the property in-place
  std::swap(prop, new_prop);

  return rts;
}

void prop_in_trans(TransitionSystem & ts, const Term & prop)
{
  // NOTE: CRUCIAL that we pass false here
  // cannot add to init or the next states
  // passing false prevents that
  ts.add_constraint(prop, false);
}

TransitionSystem promote_inputvars(const TransitionSystem & ts,
                                   const UnorderedTermSet & ivs_to_promote)
{
  SmtSolver solver = ts.solver();
  TransitionSystem new_ts = create_fresh_ts(ts.is_functional(), solver);
  assert(new_ts.is_functional() == ts.is_functional());

  // copy over all state variables
  for (const auto & sv : ts.statevars()) {
    new_ts.add_statevar(sv, ts.next(sv));
  }

  // copy over inputs but make them statevars
  for (const auto & iv : ts.inputvars()) {
    new_ts.add_inputvar(iv);
    if (ivs_to_promote.find(iv) != ivs_to_promote.end()) {
      new_ts.promote_inputvar(iv);
    }
  }

  // set init
  new_ts.set_init(ts.init());

  // copy over state updates
  for (const auto & elem : ts.state_updates()) {
    new_ts.assign_next(elem.first, elem.second);
  }

  // relational systems could have things added by constrain_trans;
  // this step has to be performed before re-adding constraints as new
  // constraints might be added to trans()
  if (!new_ts.is_functional()) {
    RelationalTransitionSystem & rts_view =
        static_cast<RelationalTransitionSystem &>(new_ts);
    rts_view.set_trans(ts.trans());
  }

  // need to re-evaluate all constraints that used to be over inputs
  for (const auto & elem : ts.constraints()) {
    new_ts.add_constraint(elem.first, elem.second);
  }

  return new_ts;
}

TransitionSystem promote_inputvars(const TransitionSystem & ts)
{
  TransitionSystem ret = promote_inputvars(ts, ts.inputvars());
  assert(ret.inputvars().empty());
  return ret;
}

}  // namespace pono
