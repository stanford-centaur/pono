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
    assert (r == ERROR);
    return "ERROR";
  }
}

}  // namespace pono
