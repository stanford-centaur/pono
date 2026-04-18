#pragma once
#include "core/ts.h"
#include "core/rts.h"

namespace pono{

/**
 * Strict: delays are > 0; Weak: delays are >=0.
*/
enum class TADelayStrictness
{
    Strict,
    Weak
};

std::string to_string(TADelayStrictness strictness);

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
 * To allow for executions to start with a time delay, a dummy discrete transition
 * is added where all variables are unchanged but time delays are allowed (respecting
 * the invariant of the initial states).
 *
 * Delays can be integer or real according to the type of clock variables. This is
 * automatically detected by the sort of clock variables (Real or Int). 
 * @see TimedVMTEncoder
*/
class TimedTransitionSystem : public RelationalTransitionSystem {
    public:
    TimedTransitionSystem(const smt::SmtSolver & s, 
        TADelayStrictness strictness = TADelayStrictness::Weak
        ) : 
        RelationalTransitionSystem(s),
        solver_(s),
        clockinvar_(solver_->make_term(true)),
        urgent_(solver_->make_term(false)),
        encoded_delays_(false),
        has_dummy_init_transitions_(false),
        delay_sort_(s->make_sort(smt::REAL)),
        delay_strictness_(strictness),
        invariant_(s->make_term(true))
        {
        }
    
    // TODO Add a copy constructor
    TimedTransitionSystem(const TimedTransitionSystem & other) = delete;

    static const std::string DELAY_VAR_NAME;
    static const std::string DUMMY_INIT_VAR_NAME;

    const smt::UnorderedTermSet & nonclock_vars() {
        return nonclock_vars_;
    }
    const smt::UnorderedTermSet & clock_vars() {
        return clock_vars_;
    }
    void add_clock_var(smt::Term var){
        clock_vars_.insert(var);
    }
    void add_nonclock_var(smt::Term var){
        nonclock_vars_.insert(var);
    }
    const smt::Term & clockinvar(){
        return clockinvar_;
    }
    const smt::Term & urgent(){
        return urgent_;
    }
    /**
     * Add an invariant constraint: the global invariant is the conjunction
     * of all added invariants.
     * 
     */
    void add_invar(const smt::Term & inv) {
        TransitionSystem::add_invar(inv);
        // TODO check if inv contains clock variables
        // TODO apply syntactic checks to clock invariants
        // TODO store clock invariants in clockinvar_
        // (we add inv in all cases using add_invar)
        clockinvar_ = solver_->make_term(smt::And, clockinvar_, inv);
    }
    /**
     * Add an urgency constraint: the global urgency constraint is the disjunction
     * of all added constraints, which all describe urgent states.
     */
    void add_urgent(const smt::Term & u){
        urgent_ = solver_->make_term(smt::Or, urgent(), u);
    }
    /**
     * Redefine the transition relation by adding delays
     * @pre urgent() only contains nonclock variables
     * @pre locinvar() only contains upper bounds on clocks (when put in negation normal form)
     * @pre invar() does not contain clock variables
    */
    void encode_timed_automaton_delays();

    protected:

    /** Add trivial edges on init states: trans := trans \/ X=X'/\ C=C' /\ init.
     * This is necessary to allow executions to start with a delay. 
     */
    void add_dummy_init_transitions();

    /**
     * Given trans(C, X, I, C', X'), add a fresh delay variable, and redefine the trans relation as
     * 
     * T(C, X, I, C',X') =
     *    C >= 0 
     * /\ delta >= 0 
     * /\ locinvar(C,X)
     * /\ (urgent -> delta = 0) 
     * /\ (trans(C,X,I,C'-delta,X')
     * /\ locinvar(X',C')
     * 
     * Moreover, to allow delays at the first step, we add a dummy init edge (without guard or reset).
     * @pre locinvar only contains upper bounds on clocks (when put in negation normal form)
    */
    void encode_compact_delays();

    const smt::SmtSolver & solver_;
    smt::UnorderedTermSet nonclock_vars_;
    smt::UnorderedTermSet clock_vars_;
    smt::Term clockinvar_;
    smt::Term urgent_;
    smt::Term delta_;
    smt::Sort delay_sort_;
    TADelayStrictness delay_strictness_;

    bool encoded_delays_;
    bool has_dummy_init_transitions_;

    smt::Term invariant_;    
};
}