#include "available_solvers.h"

// Always include boolector
#include "smt-switch/boolector_factory.h"

#if WITH_CVC4
#include "smt-switch/cvc4_factory.h"
#endif

#if WITH_MSAT
#include "smt-switch/msat_factory.h"
#endif

using namespace smt;

namespace cosa_tests {

const std::vector<SolverEnum> solver_enums({
  BTOR,

#if WITH_CVC4
      CVC4,
#endif

#if WITH_MSAT
      MSAT,
#endif
});

const CreateSolverFunsMap solvers({
  { BTOR, BoolectorSolverFactory::create },

#if WITH_CVC4
      { CVC4, CVC4SolverFactory::create },
#endif

#if WITH_MSAT
      { MSAT, MsatSolverFactory::create },
#endif
});

CreateSolverFunsMap available_solvers() { return solvers; }

std::vector<SolverEnum> available_solver_enums() { return solver_enums; }

std::ostream & operator<<(std::ostream & o, SolverEnum e)
{
  switch (e) {
    case BTOR: o << "BTOR"; break;
    case CVC4: o << "CVC4"; break;
    case MSAT: o << "MSAT"; break;
    default:
      // should print the integer representation
      throw NotImplementedException("Unknown SolverEnum: " + std::to_string(e));
      break;
  }

  return o;
}

}  // namespace cosa_tests
