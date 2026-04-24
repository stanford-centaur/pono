#pragma once
#include "core/rts.h"
#include "core/ts.h"
#include "smt-switch/utils.h"
#include "utils/exceptions.h"

namespace pono {

/**
 * Strict: delays are > 0; Weak: delays are >=0.
 */
enum class TimedAutomatonDelayStrictness
{
  Strict,
  Weak
};

std::string to_string(TimedAutomatonDelayStrictness strictness);

/**
 * Relational transition relation encoding timed automata semantics
 * where each transition is a discrete-step followed by an arbitrary time delay.
 *
 * The encoding is built in two steps. First, the front-end (e.g. VMT or SMV)
 * specifies the transition relation of the discrete steps and invariants to be
 * respected on clock variables.
 * Second, the constructor below builds (via encode_timed_automaton_delays) the
 * full transition relation by defining a single transition as a discrete step
 * followed by a timed delay along which state variables except for clocks do
 * not change, and clocks all grow by some amount, while respecting location
 * invariants and urgency constraints.
 *
 * To allow for executions to start with a time delay, a dummy discrete
 * transition is added where all variables are unchanged but time delays are
 * allowed (respecting the invariant of the initial states).
 *
 * Delays can be integer or real according to the type of clock variables. This
 * is automatically detected by the sort of clock variables (Real or Int).
 * @see TimedVMTEncoder
 */
class TimedTransitionSystem : public RelationalTransitionSystem
{
 public:
  TimedTransitionSystem() : RelationalTransitionSystem() {}
  TimedTransitionSystem(const smt::SmtSolver & s,
                        TimedAutomatonDelayStrictness strictness =
                            TimedAutomatonDelayStrictness::Weak)
      : RelationalTransitionSystem(s),
        clockinvar_(solver_->make_term(true)),
        urgent_(solver_->make_term(false)),
        encoded_delays_(false),
        delay_sort_(s->make_sort(smt::REAL)),
        delay_strictness_(strictness)
  {
  }
  TimedTransitionSystem(const TimedTransitionSystem & other,
                        smt::TermTranslator & tt);

  TimedTransitionSystem & operator=(TimedTransitionSystem other);
  friend void swap(TimedTransitionSystem & ts1, TimedTransitionSystem & ts2);

  static const std::string DELAY_VAR_NAME;

  const smt::UnorderedTermSet & nonclock_vars() { return nonclock_vars_; }
  const smt::UnorderedTermSet & clock_vars() { return clock_vars_; }
  void add_clock_var(smt::Term var) { clock_vars_.insert(var); }
  void add_nonclock_var(smt::Term var) { nonclock_vars_.insert(var); }
  const smt::Term & clockinvar() { return clockinvar_; }
  const smt::Term & urgent() { return urgent_; }
  /**
   * Add an invariant constraint: the global invariant is the conjunction
   * of all added invariants.
   *
   */
  void add_invar(const smt::Term & inv)
  {
    TransitionSystem::add_invar(inv);
    if (contains_clocks(inv)) {
      clockinvar_ = solver_->make_term(smt::And, clockinvar_, inv);
    }
  }
  /**
   * Add an urgency constraint: the global urgency constraint is the disjunction
   * of all added constraints, which all describe urgent states.
   * @pre constraint does not contain clocks.
   */
  void add_urgent(const smt::Term & u)
  {
    if (contains_clocks(u)) {
      throw PonoException("Urgency constraints cannot contain clock variables"
                          + u->to_string());
    }
    urgent_ = solver_->make_term(smt::Or, urgent(), u);
  }
  /**
   * Redefine the transition relation by adding delays.
   * To be called once after discrete transitions have been defined.
   */
  void encode_timed_automaton_delays();

 protected:
  /** Add trivial edges on init states: trans := trans \/ X=X'/\ C=C' /\ init.
   * This is necessary to allow executions to start with a delay.
   */
  void add_dummy_init_transitions();

  /**
   * Given trans(C, X, I, C', X'), add a fresh delay variable, and redefine the
   * trans relation as
   *
   * T(C, X, I, C',X') =
   *    C >= 0
   * /\ delta >= 0
   * /\ locinvar(C,X)
   * /\ (urgent -> delta = 0)
   * /\ (trans(C,X,I,C'-delta,X')
   * /\ locinvar(X',C')
   *
   * Moreover, to allow delays at the first step, we add a dummy init edge
   * (without guard or reset).
   * @pre locinvar only contains upper bounds on clocks (when put in negation
   * normal form)
   */
  void encode_compact_delays();

  /**
   * Whether the term contains clock variables.
   */
  bool contains_clocks(const smt::Term & term) const;

  smt::UnorderedTermSet nonclock_vars_;
  smt::UnorderedTermSet clock_vars_;
  smt::Term clockinvar_;
  smt::Term urgent_;
  smt::Term delta_;
  smt::Sort delay_sort_;
  TimedAutomatonDelayStrictness delay_strictness_;

  bool encoded_delays_;
};
}  // namespace pono