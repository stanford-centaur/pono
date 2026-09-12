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

/** The outcome of one run of the pono executable. */
struct PonoRun
{
  string output;   ///< standard output and standard error, interleaved
  bool succeeded;  ///< whether the process exited with a status of zero
};

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

  /** Runs pono and returns what it printed, with the standard error stream
   *  redirected into the output, since results go to standard output but the
   *  diagnostics that accompany them go to standard error.
   *  @param args the command line arguments to pass
   *  @return the run's output and whether it exited cleanly
   */
  static PonoRun run_pono(const vector<string> & args)
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
    // Zero is the only wait status standing for a clean exit, so comparing
    // against it covers an error result and a fatal signal alike.
    return { output, pclose(pipe) == 0 };
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

  /** Asserts that pono rejected a command line and said why.
   *
   *  The message reaches the output in either build, by different routes:
   *  pono.cpp's top-level try/catch prints it and an error result, but that
   *  handler is guarded by NDEBUG so a debugger or sanitizer sees the
   *  exception, so in a debug build it escapes and the terminate handler
   *  prints the same message as the process aborts. Only the printed result
   *  is therefore specific to one build.
   */
  static void expect_rejected(const PonoRun & run, const string & message)
  {
    EXPECT_FALSE(run.succeeded) << "pono exited cleanly:\n" << run.output;
    EXPECT_TRUE(contains(run.output, message));
  }
};

// The justice condition alone is violated, so the property only holds if the
// fair line is honored. k-induction is needed because bmc, the default
// engine, cannot prove a property.
TEST_F(JusticeCliUnitTests, JusticeHonorsFairnessConstraints)
{
  const PonoRun run = run_pono({ "--justice",
                                 "--engine",
                                 "ind",
                                 input_path("btor2/fair_frozen_mode.btor2") });
  EXPECT_TRUE(contains(run.output, "unsat\nj0"));
}

// Fair lines do not constrain safety checking, so checking a bad property has
// to report that they are being ignored.
TEST_F(JusticeCliUnitTests, FairnessIgnoredWithoutJustice)
{
  const PonoRun run =
      run_pono({ input_path("btor2/fair_with_bad_property.btor2") });
  EXPECT_TRUE(contains(run.output, "ignoring 1 fair line"));
}

// k-liveness counts one condition, so it has to reject a fairness constraint
// rather than silently check the justice condition alone.
TEST_F(JusticeCliUnitTests, KLivenessRejectsFairness)
{
  // Matched on enough of the message to tell it apart from a path, which can
  // itself contain the word "fairness".
  expect_rejected(run_pono({ "--justice",
                             "--justice-translator",
                             "klive",
                             input_path("btor2/fair_frozen_mode.btor2") }),
                  "including fairness constraints");
}

// Errors label the property by the kind that was requested, so a failure
// under --justice reports j, not b.
TEST_F(JusticeCliUnitTests, ErrorLabelsPropertyByRequestedKind)
{
  const PonoRun run = run_pono({ "--justice",
                                 "--prop",
                                 "5",
                                 input_path("btor2/fair_frozen_mode.btor2") });
  expect_rejected(run, "Property index 5");
#ifdef NDEBUG
  // Only the build that catches the exception reaches the result printing.
  EXPECT_TRUE(contains(run.output, "error\nj5"));
#endif
}

// The SMV and VMT encoders parse no liveness properties, so --justice has to
// be refused instead of silently checking an invariant property.
TEST_F(JusticeCliUnitTests, JusticeRejectedForSmv)
{
  expect_rejected(run_pono({ "--justice", input_path("smv/counter.smv") }),
                  "--justice is not supported for smv");
}

}  // namespace pono_tests
