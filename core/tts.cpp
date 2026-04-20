#include "core/tts.h"

#include <cassert>

#include "smt-switch/term.h"
#include "utils/logger.h"
using namespace smt;
namespace pono {
const std::string TimedTransitionSystem::DELAY_VAR_NAME = "_pono_time_delay_";

std::string to_string(TimedAutomatonDelayStrictness strictness)
{
  switch (strictness) {
    case TimedAutomatonDelayStrictness::Strict: return "Strict";
    case TimedAutomatonDelayStrictness::Weak: return "Weak";
  }
  return "Unknown";
}

void TimedTransitionSystem::encode_timed_automaton_delays()
{
  if (encoded_delays_) return;
  this->encoded_delays_ = true;
  encode_compact_delays();
}

void TimedTransitionSystem::add_dummy_init_transitions()
{
  assert(!this->has_dummy_init_transitions_);
  this->has_dummy_init_transitions_ = true;

  smt::Term dummy_transition = solver_->make_term(true);
  for (auto v : statevars_) {
    smt::Term vunchanged = solver_->make_term(Equal, v, next(v));
    dummy_transition = solver_->make_term(And, dummy_transition, vunchanged);
  }
  dummy_transition = solver_->make_term(And, dummy_transition, init());
  logger.log(4, "Dummy transitions: {}", dummy_transition);
  // dummy_transition = X'=X /\ C'=C /\ init
  // trans := trans \/ dummy_transition
  set_trans(solver_->make_term(Or, trans(), dummy_transition));
}

bool TimedTransitionSystem::contains_clocks(const smt::Term & term) const
{
  smt::UnorderedTermSet vars;
  smt::UnorderedTermSet clock_vars;
  get_free_symbolic_consts(term, vars);
  for (auto c : clock_vars_) {
    if (vars.find(c) != vars.end()) {
      return true;
    }
  }
  return false;
}

bool TimedTransitionSystem::check_clock_invariant(const smt::Term & inv) const
{
  smt::Op op = inv->get_op();
  if (op.is_null()) {
    // inv contains at least one clock variable, so cannot be a Boolean
    // constraint
    return false;
  }
  TermVec children(inv->begin(), inv->end());
  if (op.prim_op == Implies) {
    return !contains_clocks(children[0]) && is_clock_guard(children[1]);
  }
  return is_clock_guard(inv);
}

bool TimedTransitionSystem::is_clock_guard(const smt::Term & term) const
{
  if (term == solver_->make_term(true) || term == solver_->make_term(false))
    return true;
  smt::Op op = term->get_op();
  TermVec children(term->begin(), term->end());
  switch (op.prim_op) {
    case And:
      for (auto ch : children) {
        if (!is_clock_guard(ch)) {
          return false;
        }
      }
    default: return is_clock_predicate(term);
  }
}

bool TimedTransitionSystem::is_clock_predicate(const smt::Term & term) const
{
  smt::Op op = term->get_op();
  if (op.is_null()) {
    return false;
  }
  TermVec children(term->begin(), term->end());
  switch (op.prim_op) {
    case Le:
    case Lt:
    case Ge:
    case Gt:
    case Equal:
      for (auto ch : children) {
        if (!is_clock_expression(ch)) {
          return false;
        }
      }
      break;
  }
  return true;
}

bool TimedTransitionSystem::is_clock_expression(const smt::Term & term) const
{
  smt::Op op = term->get_op();
  if (op.is_null()) {
    if (term->is_value() || clock_vars_.find(term) != clock_vars_.end())
      return true;
  }
  TermVec children(term->begin(), term->end());
  if (op.prim_op == To_Real || op.prim_op == To_Int) {
    assert(children.size() == 1);
    if (children[0]->is_value()) return true;
  }
  if (op.prim_op == Minus || op.prim_op == Plus) {
    assert(children.size() == 2);
    // x - y is OK x + y is not
    if (clock_vars_.find(children[0]) != clock_vars_.end()
        && clock_vars_.find(children[1]) != clock_vars_.end()) {
      return op.prim_op == Minus;
    }
    // x +/- k
    if (children[0]->is_value()
        && clock_vars_.find(children[1]) != clock_vars_.end())
      return true;
    // k +/- x
    if (children[1]->is_value()
        && clock_vars_.find(children[0]) != clock_vars_.end())
      return true;
    // k +/- k'
    if (children[0]->is_value() && children[1]->is_value()) return true;
  }
  return false;
}

void TimedTransitionSystem::encode_compact_delays()
{
  add_dummy_init_transitions();
  smt::Sort realsort = solver_->make_sort(smt::REAL);
  smt::Sort intsort = solver_->make_sort(smt::INT);

  // determine the sort of clocks (and thus that of delays): int or real
  // (the default sort is considered if there are no clocks)
  if (clock_vars_.size() > 0) {
    smt::Term clock = *clock_vars_.begin();
    this->delay_sort_ = clock->get_sort();
    if (this->delay_sort_ != realsort && this->delay_sort_ != intsort) {
      throw PonoException("Clocks must be of sort Int or Real: "
                          + clock->to_string());
    }
    for (auto c : clock_vars_) {
      if (c->get_sort() != this->delay_sort_) {
        throw PonoException("All clocks must have the same sort: clock "
                            + c->to_string() + " is not of sort "
                            + c->get_sort()->to_string());
      }
    }
  }

  delta_ = TransitionSystem::make_inputvar(DELAY_VAR_NAME, this->delay_sort_);
  logger.log(1, "Sort of timed automata clocks: {}", this->delay_sort_);
  logger.log(
      1, "Delays are: {}", to_string(TimedAutomatonDelayStrictness::Strict));
  logger.log(2, "Listing timed automata clocks:");
  for (auto c : clock_vars_) {
    logger.log(2, "\t{}", c);
  }

  smt::Term zero, one;
  if (delay_sort_->get_sort_kind() == smt::REAL) {
    zero = solver_->make_term(To_Real,
                              solver_->make_term("0", solver_->make_sort(INT)));
    one = solver_->make_term(To_Real,
                             solver_->make_term("1", solver_->make_sort(INT)));
  } else {
    zero = solver_->make_term("0", solver_->make_sort(INT));
    one = solver_->make_term("1", solver_->make_sort(INT));
  }
  smt::Term clocks_nonnegative = solver_->make_term(true);
  smt::Term clocks_are_zero = solver_->make_term(true);
  for (auto c : clock_vars_) {
    clocks_are_zero = solver_->make_term(
        And, clocks_are_zero, solver_->make_term(Equal, c, zero));
    clocks_nonnegative = solver_->make_term(
        And, clocks_nonnegative, solver_->make_term(Ge, c, zero));
  }

  /*
   * newtrans(C, X, I, C',X') =
   *  /\ C >= 0 /\ C' >= 0
   *  /\ clockinvar(X,C)
   *  /\ delta >= 0 or > 0
   *  /\ (urgent(X') -> delta = 0)
   *  /\ (trans(C,X,I,C'-delta,X')
   *  /\ clockinvar(X',C')
   */

  // /\ C >= 0 /\ C' >= 0
  add_invar(clocks_nonnegative);

  smt::Term delta_nonnegative;
  if (this->delay_strictness_ == TimedAutomatonDelayStrictness::Weak)
    delta_nonnegative = solver_->make_term(Le, zero, delta_);
  else
    delta_nonnegative = solver_->make_term(Lt, zero, delta_);
  logger.log(4, "TA nonnegative: {}", delta_nonnegative);

  // /\ delta >= 0
  smt::Term new_trans = delta_nonnegative;

  // /\ (urgent(X') -> delta = 0)
  smt::Term delta0ifurgent = solver_->make_term(
      Implies, next(urgent()), solver_->make_term(Equal, delta_, zero));
  new_trans = solver_->make_term(And, new_trans, delta0ifurgent);
  logger.log(4, "TA delta0ifurgent: {}", delta0ifurgent);

  smt::UnorderedTermMap Cp2Cp_minus_delta;
  for (auto c : clock_vars_) {
    Cp2Cp_minus_delta[next(c)] = solver_->make_term(Minus, next(c), delta_);
  }
  logger.log(4,
             "TA Cp2Cp_minus_delta: {}",
             solver_->substitute(trans(), Cp2Cp_minus_delta));

  //  /\ (trans(C, X, I, C'-delta, X')
  new_trans = solver_->make_term(
      And, new_trans, solver_->substitute(trans(), Cp2Cp_minus_delta));

  // /\ clockinvar(X', C')
  new_trans = solver_->make_term(And, new_trans, next(this->clockinvar()));

  set_trans(new_trans);

  // init = init /\ C = 0
  set_init(solver_->make_term(And, init(), clocks_are_zero));

  logger.log(4, "TA Init: {}", init());
  logger.log(4, "TA clockinvar: {}", clockinvar());
  logger.log(4, "TA urgent: {}", urgent());
}
}  // namespace pono