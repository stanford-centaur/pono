#include "core/tts.h"

#include <cassert>

#include "smt-switch/term.h"
#include "utils/logger.h"
using namespace smt;
namespace pono {
const std::string TimedTransitionSystem::DELAY_VAR_NAME = "_pono_time_delay_";

std::string to_string(TADelayStrictness strictness)
{
  switch (strictness) {
    case TADelayStrictness::Strict: return "Strict";
    case TADelayStrictness::Weak: return "Weak";
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
  // dummy_transition = X'=X/\ C'=C /\ init
  // trans := trans \/ dummy_transition
  set_trans(solver_->make_term(Or, trans(), dummy_transition));
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
  logger.log(1, "Delays are: {}", to_string(TADelayStrictness::Strict));
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
   *  /\ delta >= 0 or > 0
   *  /\ (urgent(X') -> delta = 0)
   *  /\ (trans(C,X,I,C'-delta,X')
   *  /\ clockinvar(X',C')
   */

  // /\ C >= 0 /\ C' >= 0
  add_invar(clocks_nonnegative);

  smt::Term delta_nonnegative;
  if (this->delay_strictness_ == TADelayStrictness::Weak)
    delta_nonnegative = solver_->make_term(Le, zero, delta_);
  else 
    delta_nonnegative = solver_->make_term(Lt, zero, delta_);
  logger.log(4, "TA nonnegative: {}", delta_nonnegative);

  // /\ delta >= 0
  // new_trans = solver_->make_term(And, new_trans, delta_nonnegative);
  smt::Term new_trans = delta_nonnegative;

  // /\ (urgent(X') -> delta = 0)
  smt::Term delta0ifurgent = 
    solver_->make_term(
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

  // /\ locinvar(X', C')
  // new_trans = solver_->make_term(And, new_trans, next(locinvar()));
  new_trans = solver_->make_term(And, new_trans, next(this->clockinvar()));

  set_trans(new_trans);

  // init = init /\ C = 0
  set_init(solver_->make_term(And, init(), clocks_are_zero));

  logger.log(4, "TA Init: {}", init());
  logger.log(4, "TA clockinvar: {}", clockinvar());
  logger.log(4, "TA urgent: {}", urgent());
}
}  // namespace pono