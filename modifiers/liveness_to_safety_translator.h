#include <string>

#include "core/ts.h"
#include "smt-switch/smt.h"

namespace pono {
class LivenessToSafetyTranslator
{
 public:
  LivenessToSafetyTranslator(std::string var_prefix = "__pono_generated__");

  smt::Term translate(TransitionSystem & ts, smt::TermVec justice_conditions);

 private:
  std::string prefix_;
};
}  // namespace pono
