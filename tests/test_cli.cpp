// Tests that drive the pono executable itself, covering the wiring in
// pono.cpp -- option handling, property selection and the printed result --
// which the library-level tests cannot reach. A case is named after the input
// format it is about, and the ones holding for any format come first.

#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "utils/str_util.h"

using namespace std;

namespace pono_tests {

/** The outcome of one run of the pono executable. */
struct PonoRun
{
  string output;   ///< standard output and standard error, interleaved
  bool succeeded;  ///< whether the process exited with a status of zero
};

class CliUnitTests : public ::testing::Test
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

  static ::testing::AssertionResult omits(const string & output,
                                          const string & unexpected)
  {
    if (output.find(unexpected) == string::npos) {
      return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure()
           << "expected not to find \"" << unexpected << "\" in output:\n"
           << output;
  }

  static ::testing::AssertionResult starts_with(const string & output,
                                                const string & expected)
  {
    if (output.rfind(expected, 0) == 0) {
      return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure()
           << "expected output to start with \"" << expected << "\":\n"
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

  /** Runs pono so that it dumps a trace to a file, and returns what it
   *  wrote.
   *  @param args the command line arguments, without the dumping ones
   *  @param option the option naming the file pono should write
   *  @param extension the suffix to give that file
   *  @return the contents of the dumped file
   */
  static string dumped(vector<string> args,
                       const string & option,
                       const string & extension)
  {
    const string path = testing::TempDir() + "pono_test" + extension;
    args.insert(args.begin(), { "--witness", option, path });
    const PonoRun run = run_pono(args);
    EXPECT_TRUE(contains(run.output, "sat")) << "nothing was refuted";

    ifstream file(path);
    EXPECT_TRUE(file.is_open()) << "nothing was written to " << path;
    const string contents((istreambuf_iterator<char>(file)),
                          istreambuf_iterator<char>());
    file.close();
    remove(path.c_str());
    return contents;
  }

  /** Runs pono so that it writes a witness, and returns what it wrote.
   *  @param args the command line arguments, without the witness ones
   *  @return the contents of the witness file
   */
  static string witness(vector<string> args)
  {
    return dumped(args, "--dump-btor2-witness", ".btor2wit");
  }

  /** Runs pono so that it writes a waveform, and returns what it wrote.
   *  @param args the command line arguments, without the waveform ones
   *  @return the contents of the dumped file
   */
  static string waveform(vector<string> args)
  {
    return dumped(args, "--vcd", ".vcd");
  }

  /** Runs pono so that it writes a witness, and returns how the witness
   *  names the property it refutes.
   *  @param args the command line arguments, without the witness ones
   *  @return the second line of the witness file, after the sat it opens with
   */
  static string witness_property(vector<string> args)
  {
    istringstream lines(witness(args));
    string result;
    string property;
    getline(lines, result);
    getline(lines, property);
    EXPECT_EQ(result, "sat");
    return property;
  }
};

// The declared version comes first, ahead of the commit it was built from.
TEST_F(CliUnitTests, VersionIsReportedFirst)
{
  // PONO_VERSION is a macro set using CMake PONO_RELEASE_VERSION.
  EXPECT_TRUE(starts_with(run_pono({ "--version" }).output, PONO_VERSION));
}

// Options naming a variable take the name the file gave it, not the one the
// encoder generates from the line that declares it.
TEST_F(CliUnitTests, Btor2ResetNamesAVariableFromTheFile)
{
  const PonoRun run =
      run_pono({ "--reset", "s", input_path("btor2/input_in_bad.btor2") });
  EXPECT_EQ(run.output, "unknown\nb0\n");
}

// A name no variable has still has to be reported rather than passed over.
TEST_F(CliUnitTests, Btor2ResetRejectsAnUnknownName)
{
  const PonoRun run = run_pono(
      { "--reset", "nosuchsignal", input_path("btor2/input_in_bad.btor2") });
  expect_rejected(run, "Could not find term named: nosuchsignal");
}

// The witness file names the property the same way the result does.
TEST_F(CliUnitTests, Btor2WitnessNamesTheProperty)
{
  EXPECT_EQ(witness_property({ input_path("btor2/input_in_bad.btor2") }), "b0");
  EXPECT_EQ(
      witness_property({ "--justice", input_path("btor2/justice_only.btor2") }),
      "j0");
}

// Pono adds variables to the system it checks for its own bookkeeping: the
// property monitor and pseudo initial state here, the lasso and the counters
// elsewhere. None of them belong to the design, so no trace reports them.
TEST_F(CliUnitTests, Btor2WaveformLeavesOutGeneratedSignals)
{
  const string vcd = waveform({ "--pseudo-init-prop",
                                "-k",
                                "3",
                                input_path("btor2/input_in_bad.btor2") });
  EXPECT_TRUE(omits(vcd, pono::generated_prefix));
  EXPECT_TRUE(contains(vcd, " s $end")) << "the design's own state is gone";
}

// Both translators add state, and neither reaches the waveform.
TEST_F(CliUnitTests, Btor2JusticeWaveformLeavesOutGeneratedSignals)
{
  EXPECT_TRUE(
      omits(waveform({ "--justice", input_path("btor2/justice_only.btor2") }),
            pono::generated_prefix));
  EXPECT_TRUE(omits(waveform({ "--justice",
                               "--justice-translator",
                               "klive",
                               input_path("btor2/justice_only.btor2") }),
                    pono::generated_prefix));
}

TEST_F(CliUnitTests, Btor2WitnessLeavesOutGeneratedSignals)
{
  EXPECT_TRUE(
      omits(witness({ "--justice", input_path("btor2/justice_only.btor2") }),
            pono::generated_prefix));
}

// The plain text trace, which the formats without a witness form fall back
// on, reports the design's variables only. That drops what pono generated
// along with the next-state twin of each state, whose value is just what the
// following step already reports.
TEST_F(CliUnitTests, SmvTraceLeavesOutGeneratedSignals)
{
  const PonoRun run = run_pono(
      { "--witness", "-k", "5", input_path("smv/counter_bitvector.smv") });
  EXPECT_TRUE(contains(run.output, "counter : "));
  EXPECT_TRUE(omits(run.output, pono::generated_prefix));
}

// A design that names a variable the way pono names its own would have that
// variable left out of every trace, so the name is refused instead.
TEST_F(CliUnitTests, Btor2RejectsADesignNamedLikeGenerated)
{
  expect_rejected(run_pono({ input_path("btor2/generated_name_clash.btor2") }),
                  "is named like one pono generates for itself");
}

// A justice counterexample is a lasso, and its prefix dumps as a waveform
// like any other trace. Both translators reach the same printer.
TEST_F(CliUnitTests, Btor2JusticeWaveformHoldsTheDesign)
{
  EXPECT_TRUE(contains(
      waveform({ "--justice", input_path("btor2/justice_only.btor2") }),
      "mode"));
}

TEST_F(CliUnitTests, Btor2JusticeWaveformHoldsTheDesignWithKLiveness)
{
  EXPECT_TRUE(contains(waveform({ "--justice",
                                  "--justice-translator",
                                  "klive",
                                  input_path("btor2/justice_only.btor2") }),
                       "mode"));
}

// A justice property with no fairness constraint: the lasso at mode=1
// satisfies the condition infinitely often, so it is a counterexample.
TEST_F(CliUnitTests, Btor2JusticeOnly)
{
  const PonoRun run =
      run_pono({ "--justice", input_path("btor2/justice_only.btor2") });
  EXPECT_EQ(run.output, "sat\nj0\n");
}

// The same property under the other translator, which counts observations of
// the single condition rather than translating it to a safety property.
TEST_F(CliUnitTests, Btor2JusticeOnlyWithKLiveness)
{
  const PonoRun run = run_pono({ "--justice",
                                 "--justice-translator",
                                 "klive",
                                 input_path("btor2/justice_only.btor2") });
  EXPECT_EQ(run.output, "sat\nj0\n");
}

// A justice property that holds, which is the outcome the cases above cannot
// reach: no lasso satisfies the condition infinitely often.
TEST_F(CliUnitTests, Btor2JusticeHolds)
{
  const PonoRun run = run_pono({ "--justice",
                                 "--engine",
                                 "ind",
                                 "--bound",
                                 "20",
                                 input_path("btor2/justice_holds.btor2") });
  EXPECT_EQ(run.output, "unsat\nj0\n");
}

// Each justice line is a separate property, so --prop picks between them. The
// two here disagree, which is what makes the selection visible.
TEST_F(CliUnitTests, Btor2SecondJusticeProperty)
{
  const PonoRun run =
      run_pono({ "--justice",
                 "--prop",
                 "1",
                 "--engine",
                 "ind",
                 "--bound",
                 "20",
                 input_path("btor2/justice_two_properties.btor2") });
  EXPECT_EQ(run.output, "unsat\nj1\n");
}

// Reducing to the state the justice condition depends on leaves the verdict
// alone.
TEST_F(CliUnitTests, Btor2JusticeWithConeOfInfluence)
{
  const PonoRun run =
      run_pono({ "--justice",
                 "--static-coi",
                 input_path("btor2/justice_unrelated_state.btor2") });
  EXPECT_EQ(run.output, "sat\nj0\n");
}

// k-liveness takes its own cone of influence rather than going through the
// one the driver takes, so the reduction reaches it by a separate route.
TEST_F(CliUnitTests, Btor2JusticeWithConeOfInfluenceAndKLiveness)
{
  const PonoRun run = run_pono({ "--justice",
                                 "--justice-translator",
                                 "klive",
                                 "--static-coi",
                                 input_path("btor2/justice_only.btor2") });
  EXPECT_EQ(run.output, "sat\nj0\n");
}

// The justice condition alone is violated, so the property only holds if the
// fair line is honored. k-induction is needed because bmc, the default
// engine, cannot prove a property.
TEST_F(CliUnitTests, Btor2JusticeHonorsFairnessConstraints)
{
  const PonoRun run = run_pono({ "--justice",
                                 "--engine",
                                 "ind",
                                 input_path("btor2/fair_frozen_mode.btor2") });
  EXPECT_TRUE(contains(run.output, "unsat\nj0"));
}

// Btor2 fairness constraints apply to justice properties only, so a fair line
// has no bearing on checking a bad property. The whole output is matched
// rather than searched, because "sat" occurs in "unsat" too.
TEST_F(CliUnitTests, Btor2FairnessDoesNotAffectSafetyChecking)
{
  const PonoRun run =
      run_pono({ input_path("btor2/fair_with_bad_property.btor2") });
  EXPECT_EQ(run.output, "sat\nb0\n");
}

// The cone of influence is taken over the fairness constraints as well, so
// the state only they watch survives it. Reducing over the justice condition
// alone drops that state and the run dies translating the property.
TEST_F(CliUnitTests, Btor2JusticeWithConeOfInfluenceKeepsFairnessState)
{
  const PonoRun run =
      run_pono({ "--justice",
                 "--static-coi",
                 input_path("btor2/fair_independent_state.btor2") });
  EXPECT_EQ(run.output, "sat\nj0\n");
}

// k-liveness counts one condition, so it has to reject a fairness constraint
// rather than silently check the justice condition alone.
TEST_F(CliUnitTests, Btor2KLivenessRejectsFairness)
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
TEST_F(CliUnitTests, Btor2ErrorLabelsPropertyByRequestedKind)
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
TEST_F(CliUnitTests, SmvJusticeIsRejected)
{
  expect_rejected(run_pono({ "--justice", input_path("smv/counter.smv") }),
                  "--justice is not supported for smv");
}
}  // namespace pono_tests
