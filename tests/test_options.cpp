#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "options/options.h"

using namespace pono;
using namespace std;

namespace pono_tests {

// The vector overload synthesizes an argv whose first entry is a dummy
// program name, because the argc/argv overload skips argv[0]. Asserting on
// the first option pins that contract: were the count off by one, the parser
// would take "--bound" for the program name and reject "5" as a non-option.
TEST(OptionsUnitTests, ParseOptionsVector)
{
  PonoOptions opts;
  vector<string> args({ "--bound", "5", "--engine", "bmc" });
  EXPECT_EQ(opts.parse_and_set_options(args, false), UNKNOWN);
  EXPECT_EQ(opts.bound_, 5ul);
  EXPECT_EQ(opts.engine_, BMC);
}

}  // namespace pono_tests
