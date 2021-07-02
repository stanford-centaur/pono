#include "proverresult.h"
#include <cassert>

namespace pono {

std::string to_string(ProverResult r)
{
  if (r == TRUE) {
    return "TRUE";
  } else if (r == FALSE) {
    return "FALSE";
  } else if (r == UNKNOWN) {
    return "UNKNOWN";
  } else {
    assert(r == ERROR);
    return "ERROR";
  }
}

std::ostream & operator<<(std::ostream & o, ProverResult r)
{
  o << to_string(r);
  return o;
}

}  // namespace pono
