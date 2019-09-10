#include "prop.h"
#include "exceptions.h"

using namespace smt;

namespace cosa
{

  Property::Property(const RelationalTransitionSystem &ts, Term p):
    ts_(ts),
    prop_(p)
  {
  }

  Property::~Property()
  {
  }

} // namespace cosa
