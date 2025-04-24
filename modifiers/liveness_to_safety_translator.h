/**
 * \file modifiers/liveness_to_safety_translator.h
 * \brief Liveness to safety translation
 *
 * Based on: Armin Biere, Cyrille Artho, and Viktor Schuppan,
 *           “Liveness Checking as Safety Checking,”
 *           Electronic Notes in Theoretical Computer Science,
 *           Volume 66, Issue 2, pages 160–177, 2002.
 */

#include <string>

#include "core/ts.h"
#include "smt-switch/smt.h"

namespace pono {
class LivenessToSafetyTranslator
{
 public:
  /**
   * \brief Initialize liveness-to-safety translator.
   *
   * \param var_prefix text that will be prepended to the names of all
   *                   generated variables; can be useful to avoid name
   *                   collisions
   */
  LivenessToSafetyTranslator(std::string var_prefix = "__pono_generated__");

  /**
   * \brief Perform translation of liveness property to safety.
   *
   * The given justice properties j(i) are equivalent to the LTL property given
   * in the AIGER spec: GFj(i) must hold for all i on the witness trace.
   *
   * This will add a new oracle input "save" that marks a state which is the
   * first in the loop part of a lasso-shaped trace. A new input "saved" is
   * added to mark whether the input has been set to true in the past.
   *
   * For each existing state variable, a new "loop" copy is created to store
   * the value of each variable at the state where "save" was (first) true.
   * For each jusice condition, a boolean state variable is added that stores
   * whether it has held at any point following entering the loop.
   *
   *
   * \param ts transition system to perform the translation on
   * \param justice_conditions set of justice conditions that must hold
   *                           infinitely often in the witness trace
   * \retval safety_property becomes false if the current state is the same as
   *         a previous state (the one where "save" was first true) and every
   *         justice condition has been satisfied at least once in the loop
   */
  smt::Term translate(TransitionSystem & ts, smt::TermVec justice_conditions);

 private:
  std::string prefix_;
};
}  // namespace pono
