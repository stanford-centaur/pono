#include "engines/liveness_to_safety.h"

#include "core/prop.h"
#include "options/options.h"

namespace pono {

template <class SafetyProver_T>
LivenessToSafety<SafetyProver_T>::LivenessToSafety(
    const LivenessProperty & property,
    const TransitionSystem & ts,
    const smt::SmtSolver & solver,
    PonoOptions options,
    std::string var_prefix)
    : LivenessProver(property, ts, solver, options), prefix_(var_prefix)
{
}

}  // namespace pono
