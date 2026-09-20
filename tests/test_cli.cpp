// Tests that drive the pono executable itself, for behavior that belongs to
// no particular input format.

#include <string>

#include "cli_test_fixture.h"
#include "gtest/gtest.h"

using namespace std;

namespace pono_tests {

class GeneralCliUnitTests : public CliUnitTests
{
};

// The declared version comes first, ahead of the commit it was built from.
TEST_F(GeneralCliUnitTests, VersionIsReportedFirst)
{
  // PONO_VERSION is a macro set using CMake PONO_RELEASE_VERSION.
  EXPECT_TRUE(starts_with(run_pono({ "--version" }).output, PONO_VERSION));
}

}  // namespace pono_tests
