#include "proverresult.h"

namespace pono {

std::string to_string(ProverResult r)
{
  if (r == TRUE) {
    return "TRUE";
  } else if (r == FALSE) {
    return "FALSE";
  } else {
    return "UNKNOWN";
  }
}

}  // namespace pono
