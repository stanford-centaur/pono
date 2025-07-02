/**
 * \file engines/liveness_to_safety.h
 * \brief Liveness prover that uses direct translation to safety
 *
 * Based on: Armin Biere, Cyrille Artho, and Viktor Schuppan,
 *           “Liveness Checking as Safety Checking,”
 *           Electronic Notes in Theoretical Computer Science,
 *           Volume 66, Issue 2, pages 160–177, 2002.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "core/proverresult.h"
#include "core/ts.h"
#include "engines/prover.h"
#include "smt-switch/smt.h"

namespace pono {

/**
 * \brief This prover performs translation of the liveness property to safety.
 *
 * The given justice properties j(i) are equivalent to the LTL property given
 * in the AIGER spec: GFj(i) must hold for all i on the witness trace, i.e.,
 * each justice condition must hold infinitely often.
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
 * The created safety property becomes false if the current state is the same
 * as a previous state (the one where "save" was first true) and every justice
 * condition has been satisfied at least once in the loop*
 */
template <class SafetyProver_T>
class LivenessToSafety : public LivenessProver
{
 public:
  /**
   * \brief Initialize liveness-to-safety translation engine.
   *
   * \param var_prefix text that will be prepended to the names of all
   *                   generated variables; can be useful to avoid name
   *                   collisions
   */
  LivenessToSafety(const LivenessProperty & property,
                   const TransitionSystem & ts,
                   const smt::SmtSolver & solver,
                   PonoOptions options = {},
                   std::string var_prefix = "__pono_generated__");

  void initialize() override;

  ProverResult prove() override;

  ProverResult check_until(int k) override;

  virtual bool witness(std::vector<smt::UnorderedTermMap> & out) override;

  virtual std::size_t witness_length() const override;

 private:
  std::string prefix_;

  SafetyProver_T safety_prover_;
};  // class LivenessToSafety

}  // namespace pono
