#include "core/tts.h"

#include <cassert>
#include <functional>

#include "smt-switch/substitution_walker.h"
#include "smt-switch/term.h"
#include "smt-switch/utils.h"
#include "smt/available_solvers.h"
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
void swap(TimedTransitionSystem & ts1, TimedTransitionSystem & ts2)
{
  swap((TransitionSystem &)ts1, (TransitionSystem &)ts2);
  std::swap(ts1.clock_vars_, ts2.clock_vars_);
  std::swap(ts1.nonclock_vars_, ts2.nonclock_vars_);
  std::swap(ts1.delay_strictness_, ts2.delay_strictness_);
  std::swap(ts1.clockinvar_, ts2.clockinvar_);
  std::swap(ts1.urgent_, ts2.urgent_);
  std::swap(ts1.delta_, ts2.delta_);
  std::swap(ts1.delay_sort_, ts2.delay_sort_);
  std::swap(ts1.encoded_delays_, ts2.encoded_delays_);
}
TimedTransitionSystem & TimedTransitionSystem::operator=(
    TimedTransitionSystem other)
{
  swap(*this, other);
  return *this;
}
TimedTransitionSystem::TimedTransitionSystem(
    const TimedTransitionSystem & other_ts, smt::TermTranslator & tt)
    : RelationalTransitionSystem(other_ts, tt)
{
  std::function<Term(const Term &)> transfer;
  std::function<Term(const Term &, SortKind)> transfer_as;
  if (other_ts.solver() == tt.get_solver()) {
    // if the solvers are the same, don't need to transfer
    transfer = [](const Term & t) { return t; };
    // assume you don't need to do sort-casting for terms from the same solver
    transfer_as = [](const Term & t, SortKind sk) { return t; };
    delay_sort_ = other_ts.delay_sort_;
  } else {
    transfer = [&tt](const Term & t) { return tt.transfer_term(t); };
    transfer_as = [&tt](const Term & t, SortKind sk) {
      return tt.transfer_term(t, sk);
    };
    delay_sort_ = tt.transfer_sort(other_ts.delay_sort_);
  }

  for (const auto & v : other_ts.clock_vars_) {
    clock_vars_.insert(transfer(v));
  }
  for (const auto & v : other_ts.nonclock_vars_) {
    nonclock_vars_.insert(transfer(v));
  }
  delta_ = transfer(other_ts.delta_);
  clockinvar_ = transfer_as(other_ts.clockinvar_, BOOL);
  urgent_ = transfer_as(other_ts.urgent_, BOOL);
  delay_strictness_ = other_ts.delay_strictness_;
  encoded_delays_ = other_ts.encoded_delays_;
}

void TimedTransitionSystem::encode_timed_automaton_delays()
{
  if (encoded_delays_)
    throw PonoException(
        "The encode_timed_automaton_delays function must be called only once.");
  this->encoded_delays_ = true;
  encode_compact_delays();
}

void TimedTransitionSystem::add_dummy_init_transitions()
{
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