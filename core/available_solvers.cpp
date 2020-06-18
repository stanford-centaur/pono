#include "core/available_solvers.h"
#include "utils/exceptions.h"

#include <sstream>

// Always include boolector
#include "smt-switch/boolector_factory.h"

#if WITH_CVC4
#include "smt-switch/cvc4_factory.h"
#endif

#if WITH_MSAT
#include "smt-switch/msat_factory.h"
#endif

#if WITH_YICES2
#include "smt-switch/yices2_factory.h"
#endif

using namespace smt;

namespace cosa {

// list of regular (non-interpolator) solver enums
const std::vector<SolverEnum> solver_enums({
  BTOR, BTOR_LOGGING,

#if WITH_CVC4
      CVC4, CVC4_LOGGING,
#endif

#if WITH_MSAT
      MSAT, MSAT_LOGGING,
#endif

#if WITH_YICES2
      YICES2, YICES2_LOGGING,
#endif
});

const std::vector<SolverEnum> itp_enums({
#if WITH_MSAT
  MSAT_INTERPOLATOR
#endif
});

SmtSolver create_solver(SolverEnum se)
{
  switch (se) {
    case BTOR: {
      return BoolectorSolverFactory::create(false);
      break;
      ;
    }
    case BTOR_LOGGING: {
      return BoolectorSolverFactory::create(true);
      break;
      ;
    }
#if WITH_CVC4
    case CVC4: {
      return CVC4SolverFactory::create(false);
      break;
      ;
    }
    case CVC4_LOGGING: {
      return CVC4SolverFactory::create(true);
      break;
      ;
    }
#endif
#if WITH_MSAT
    case MSAT: {
      return MsatSolverFactory::create(false);
      break;
      ;
    }
    case MSAT_LOGGING: {
      return MsatSolverFactory::create(true);
      break;
      ;
    }
#endif
#if WITH_YICES2
    case YICES2: {
      return Yices2SolverFactory::create(false);
      break;
      ;
    }
    case YICES2_LOGGING: {
      return Yices2SolverFactory::create(true);
      break;
      ;
    }
#endif
    default: {
      std::ostringstream oss;
      oss << "Not built with solver " << se;
      throw CosaException(oss.str());
    }
  }
}

SmtSolver create_interpolator(SolverEnum se)
{
  switch (se) {
#if WITH_MSAT
    case MSAT_INTERPOLATOR: {
      return MsatSolverFactory::create_interpolating_solver();
      break;
      ;
    }
#endif
    default: {
      std::ostringstream oss;
      oss << "Not built with interpolator " << se;
      throw CosaException(oss.str());
    }
  }
}

std::vector<SolverEnum> available_solver_enums() { return solver_enums; }

std::vector<SolverEnum> available_interpolator_enums() { return itp_enums; };

}  // namespace cosa
