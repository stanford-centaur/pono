#include "available_solvers.h"

using namespace smt;

namespace pono_tests {

const std::vector<SolverEnum> solver_enums({
  BTOR, CVC4,

#if WITH_MSAT
      MSAT,
#endif
});

const CreateSolverFunsMap solvers({
  { BTOR, BoolectorSolverFactory::create }, { CVC4, CVC4SolverFactory::create },

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

}  // namespace pono_tests
