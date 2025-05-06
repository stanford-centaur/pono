#include "modifiers/liveness_to_safety_translator.h"

#include <string>
#include <utility>
#include <vector>

#include "core/ts.h"
#include "smt-switch/smt.h"

namespace pono {
LivenessToSafetyTranslator::LivenessToSafetyTranslator(std::string var_prefix)
    : prefix_(var_prefix)
{
}

smt::Term LivenessToSafetyTranslator::translate(TransitionSystem & ts,
                                                smt::TermVec justice_conditions)
{
  auto boolsort = ts.solver()->make_sort(smt::BOOL);
  auto false_val = ts.solver()->make_term(false);
  auto true_val = ts.solver()->make_term(true);
  // Store a copy of the set of pre-existing state variables.
  const smt::UnorderedTermSet orig_statevars = ts.statevars();

  // Add oracle input. When this becomes true, we "save" the current state,
  // then continue until we find the same state again.
  auto save_input = ts.make_inputvar(prefix_ + "save", boolsort);

  // Add "saved" state. This indicates whether we have already saved a state,
  // i.e., whether we are inside a loop.
  auto saved_state = ts.make_statevar(prefix_ + "saved", boolsort);
  // Starts as false.
  ts.constrain_init(ts.solver()->make_term(smt::Equal, saved_state, false_val));
  // saved' = saved \/ save
  ts.assign_next(saved_state,
                 ts.solver()->make_term(smt::Or, saved_state, save_input));

  // Add "loop" states. These keep track of the first state in the loop.
  std::vector<std::pair<smt::Term, smt::Term>> loop_states;
  for (auto statevar : orig_statevars) {
    auto loop_state = ts.make_statevar(
        statevar->to_string() + prefix_ + "_loop", statevar->get_sort());
    loop_states.push_back({ statevar, loop_state });
    // loop_i' = (save /\ !saved) ? state_i : loop_i;
    ts.assign_next(loop_state,
                   ts.solver()->make_term(
                       smt::Ite,
                       { ts.solver()->make_term(
                             smt::And,
                             save_input,
                             ts.solver()->make_term(smt::Not, saved_state)),
                         statevar,
                         loop_state }));
  }

  // Add justice states. They start as false and become true once we are
  // inside a loop and the relevant justice property is met.
  smt::TermVec justice_states;
  for (std::size_t i = 0; i < justice_conditions.size(); i++) {
    auto justice_state =
        ts.make_statevar(prefix_ + "justice_" + std::to_string(i), boolsort);
    justice_states.push_back(justice_state);
    ts.constrain_init(
        ts.solver()->make_term(smt::Equal, justice_state, false_val));
    ts.assign_next(
        justice_state,
        ts.solver()->make_term(
            smt::Or,
            justice_state,
            ts.solver()->make_term(
                smt::And,
                ts.solver()->make_term(smt::Or, save_input, saved_state),
                justice_conditions[i])));
  }

  // Construct safety property.
  // looped := saved /\ (state_0 = loop_0) /\ ... /\ (state_i = loop_i)
  auto looped = saved_state;
  for (const auto & var_pair : loop_states) {
    looped = ts.solver()->make_term(
        smt::And,
        looped,
        ts.solver()->make_term(smt::Equal, var_pair.first, var_pair.second));
  }
  // bad := looped /\ justice_0 /\ ... /\ justice_j
  auto bad = looped;
  for (auto justice_state : justice_states) {
    bad = ts.solver()->make_term(smt::And, bad, justice_state);
  }

  // safe = !bad
  return ts.solver()->make_term(smt::Not, bad);
}
}  // namespace pono
