// Tests that drive the pono executable itself, for justice properties and the
// fairness constraints that go with them. They cover the wiring in pono.cpp --
// option handling, property selection and the printed result -- which the
// library-level tests cannot reach.

#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include "gtest/gtest.h"

using namespace std;

namespace pono_tests {

class JusticeCliUnitTests : public ::testing::Test
{
 protected:
  // PONO_INPUT_DIR is a macro set using CMake to the test input directory.
  static string input_path(const string & name)
  {
    return string(PONO_INPUT_DIR) + "/" + name;
  }

  /** Wraps a word in single quotes so the shell takes it literally. */
  static string quote(const string & word)
  {
    string quoted = "'";
    for (const char c : word) {
      // A single quote has to leave the quoted run to be escaped.
      if (c == '\'') {
        quoted += "'\\''";
      } else {
        quoted += c;
      }
    }
    return quoted + "'";
  }

  /** Runs pono and returns its output, with the standard error stream
   *  redirected into it, since results go to standard output but the
   *  diagnostics that accompany them go to standard error.
   *  @param args the command line arguments to pass
   *  @return the combined output
   */
  static string run_pono(const vector<string> & args)
  {
    // PONO_BIN is a macro set using CMake to the executable's location.
    string command = quote(PONO_BIN);
    for (const string & arg : args) {
      command += " " + quote(arg);
    }
    command += " 2>&1";

    FILE * pipe = popen(command.c_str(), "r");
    if (!pipe) {
      throw runtime_error("could not run: " + command);
    }
    string output;
    char buffer[256];
    while (fgets(buffer, sizeof(buffer), pipe)) {
      output += buffer;
    }
    pclose(pipe);
    return output;
  }

  static ::testing::AssertionResult contains(const string & output,
                                             const string & expected)
  {
    if (output.find(expected) != string::npos) {
      return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure()
           << "expected to find \"" << expected << "\" in output:\n"
           << output;
  }
};

// The justice condition alone is violated, so the property only holds if the
// fair line is honored. k-induction is needed because bmc, the default
// engine, cannot prove a property.
TEST_F(JusticeCliUnitTests, JusticeHonorsFairnessConstraints)
{
  const string output =
      run_pono({ "--justice",
                 "--engine",
                 "ind",
                 input_path("btor2/fair_frozen_mode.btor2") });
  EXPECT_TRUE(contains(output, "unsat\nj0"));
}

// Fair lines do not constrain safety checking, so checking a bad property has
// to report that they are being ignored.
TEST_F(JusticeCliUnitTests, FairnessIgnoredWithoutJustice)
{
  const string output =
      run_pono({ input_path("btor2/fair_with_bad_property.btor2") });
  EXPECT_TRUE(contains(output, "ignoring 1 fair line"));
}

// k-liveness counts one condition, so it has to reject a fairness constraint
// rather than silently check the justice condition alone.
TEST_F(JusticeCliUnitTests, KLivenessRejectsFairness)
{
  const string output =
      run_pono({ "--justice",
                 "--justice-translator",
                 "klive",
                 input_path("btor2/fair_frozen_mode.btor2") });
  EXPECT_TRUE(contains(output, "fairness"));
}

// Errors label the property by the kind that was requested, so a failure
// under --justice reports j, not b.
TEST_F(JusticeCliUnitTests, ErrorLabelsPropertyByRequestedKind)
{
  const string output =
      run_pono({ "--justice",
                 "--prop",
                 "5",
                 input_path("btor2/fair_frozen_mode.btor2") });
  EXPECT_TRUE(contains(output, "error\nj5"));
}

// The SMV and VMT encoders parse no liveness properties, so --justice has to
// be refused instead of silently checking an invariant property.
TEST_F(JusticeCliUnitTests, JusticeRejectedForSmv)
{
  const string output =
      run_pono({ "--justice", input_path("smv/counter.smv") });
  EXPECT_TRUE(contains(output, "--justice is not supported for smv"));
}

}  // namespace pono_tests
