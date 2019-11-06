#include "prover.h"

#include <climits>

namespace cosa {

Prover::Prover()
{
}

Prover::~Prover()
{
}

ProverResult Prover::prove()
{
  return check_until(INT_MAX);
}
	
} // namespace cosa
