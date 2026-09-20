// Shared plumbing for tests that drive the pono executable itself, covering
// the wiring in pono.cpp -- option handling, property selection and the
// printed result -- which the library-level tests cannot reach.

#pragma once

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
};

}  // namespace pono_tests
