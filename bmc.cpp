#include "bmc.h"

using namespace smt;

namespace cosa
{

  Bmc::Bmc(const RelationalTransitionSystem &ts, const Property &p,
	   SmtSolver &solver):
    ts_(ts),
    property_(p),
    solver_(solver),
    unroller_(ts, solver)
  {
  }

  Bmc::~Bmc()
  {
  }

  bool Bmc::check_until(size_t k)
  {
    // TODO
    return false;
  }

  bool Bmc::step()
  {
    //TODO
    return false;
  }
  
} // namespace cosa
