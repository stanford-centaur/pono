// Tests that drive the pono executable itself on Btor2 input.

#include <string>

#include "cli_test_fixture.h"
#include "gtest/gtest.h"

using namespace std;

namespace pono_tests {

class Btor2CliUnitTests : public CliUnitTests
{
};

// Options naming a variable take the name the file gave it, not the one the
// encoder generates from the line that declares it.
TEST_F(Btor2CliUnitTests, ResetNamesAVariableFromTheFile)
{
  const PonoRun run =
      run_pono({ "--reset", "s", input_path("btor2/input_in_bad.btor2") });
  EXPECT_EQ(run.output, "unknown\nb0\n");
}

// A name no variable has still has to be reported rather than passed over.
TEST_F(Btor2CliUnitTests, ResetRejectsAnUnknownName)
{
  const PonoRun run = run_pono(
      { "--reset", "nosuchsignal", input_path("btor2/input_in_bad.btor2") });
  expect_rejected(run, "Could not find term named: nosuchsignal");
}

}  // namespace pono_tests
