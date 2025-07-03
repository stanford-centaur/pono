#include "core/prop.h"

#include <string>

#include "smt-switch/smt.h"

namespace pono {

AbstractProperty::AbstractProperty(const smt::SmtSolver & solver,
                                   std::string name)
    : solver_(solver), name_(name)
{
}

const smt::SmtSolver & AbstractProperty::solver() const { return solver_; }

const std::string AbstractProperty::name() const { return name_; }

SafetyProperty::SafetyProperty(const smt::SmtSolver & solver,
                               const smt::Term & property,
                               std::string name)
    : AbstractProperty(solver, name), prop_(property)
{
}

const smt::Term & SafetyProperty::prop() const { return prop_; }

LivenessProperty::LivenessProperty(const smt::SmtSolver & solver,
                                   const smt::TermVec & conditions,
                                   std::string name)
    : AbstractProperty(solver, name), prop_terms_(conditions)
{
}

const smt::TermVec & LivenessProperty::terms() const { return prop_terms_; }

}  // namespace pono
