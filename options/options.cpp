/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Florian Lonsing, Makai Mann
 ** This file is part of the pono project.
 ** Copyright (c) 2019, 2020 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief
 **
 **
 **/

#include "options/options.h"
#include <iostream>
#include <string>
#include <vector>
#include "optionparser.h"
#include "utils/exceptions.h"

using namespace std;

/************************************* Option Handling setup
 * *****************************************/
// from optionparser-1.7 examples -- example_arg.cc
enum optionIndex
{
  UNKNOWN_OPTION,
  HELP,
  ENGINE,
  BOUND,
  PROP,
  VERBOSITY,
  VCDNAME,
  NOWITNESS,
  CEGPROPHARR,
  NO_CEGP_AXIOM_RED,
  RESET,
  RESET_BND,
  CLK
};

struct Arg : public option::Arg
{
  static void printError(const char * msg1,
                         const option::Option & opt,
                         const char * msg2)
  {
    fprintf(stderr, "%s", msg1);
    fwrite(opt.name, opt.namelen, 1, stderr);
    fprintf(stderr, "%s", msg2);
  }

  static option::ArgStatus Numeric(const option::Option & option, bool msg)
  {
    char * endptr = 0;
    if (option.arg != 0 && strtol(option.arg, &endptr, 10)) {
    };
    if (endptr != option.arg && *endptr == 0) return option::ARG_OK;

    if (msg) printError("Option '", option, "' requires a numeric argument\n");
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus NonEmpty(const option::Option & option, bool msg)
  {
    if (option.arg != 0 && option.arg[0] != 0) return option::ARG_OK;

    if (msg)
      printError("Option '", option, "' requires a non-empty argument\n");
    return option::ARG_ILLEGAL;
  }
};

const option::Descriptor usage[] = {
  { UNKNOWN_OPTION,
    0,
    "",
    "",
    Arg::None,
    "USAGE: pono [options] <btor file>\n\n"
    "Options:" },
  { HELP, 0, "", "help", Arg::None, "  --help \tPrint usage and exit." },
  { ENGINE,
    0,
    "e",
    "engine",
    Arg::NonEmpty,
    "  --engine, -e <engine> \tSelect engine from [bmc, bmc-sp, ind, "
    "interp]." },
  { BOUND,
    0,
    "k",
    "bound",
    Arg::Numeric,
    "  --bound, -k \tBound to check up until (default: 10)." },
  { PROP,
    0,
    "p",
    "prop",
    Arg::Numeric,
    "  --prop, -p \tProperty index to check (default: 0)." },
  { VERBOSITY,
    0,
    "v",
    "verbosity",
    Arg::Numeric,
    "  --verbosity, -v \tVerbosity for printing to standard out." },
  { VCDNAME,
    0,
    "",
    "vcd",
    Arg::NonEmpty,
    "  --vcd \tName of Value Change Dump (VCD) if witness exists." },
  { NOWITNESS,
    0,
    "",
    "no-witness",
    Arg::None,
    "  --no-witness \tDisable printing of witness." },
  { CEGPROPHARR,
    0,
    "",
    "ceg-prophecy-arrays",
    Arg::None,
    "  --ceg-prophecy-arrays \tUse counter-example guided prophecy for "
    "arrays." },
  { NO_CEGP_AXIOM_RED,
    0,
    "",
    "no-cegp-axiom-red",
    Arg::None,
    "  --no-cegp-axiom-red \tDon't reduce axioms in CEG-Prophecy with unsat "
    "cores." },
  { RESET,
    0,
    "r",
    "reset",
    Arg::NonEmpty,
    "  --reset, -r <reset input> \tSymbol to use for reset signal (prefix with "
    "~ "
    "for negative reset)" },
  { RESET_BND,
    0,
    "s",
    "resetsteps",
    Arg::Numeric,
    "  --resetsteps, -s <integer> \tNumber of steps to apply reset for "
    "(default: 1)" },
  { CLK,
    0,
    "c",
    "clock",
    Arg::NonEmpty,
    "  --clock, -c <clock name> \tSymbol to use for clock signal (only "
    "supports "
    "starting at 0 and toggling each step)" },
  { 0, 0, 0, 0, 0, 0 }
};
/*********************************** end Option Handling setup
 * ***************************************/

namespace pono {

Engine PonoOptions::to_engine(std::string s)
{
  if (str2engine.find(s) != str2engine.end()) {
    return str2engine.at(s);
  } else {
    throw PonoException("Unrecognized engine: " + s);
  }
}

// Parse command line options given by 'argc' and 'argv' and set
// respective options in the 'pono_options' object.
// Returns 'ERROR' if there is something wrong with the given options
// or 'UNKNOWN' instead.
ProverResult PonoOptions::parse_and_set_options(int argc, char ** argv)
{
  argc -= (argc > 0);
  argv += (argc > 0);  // skip program name argv[0] if present
  option::Stats stats(usage, argc, argv);
  std::vector<option::Option> options(stats.options_max);
  std::vector<option::Option> buffer(stats.buffer_max);
  option::Parser parse(usage, argc, argv, &options[0], &buffer[0]);

  if (parse.error()) return ERROR;

  if (options[HELP] || argc == 0) {
    option::printUsage(cout, usage);
    // want to exit main at top-level
    return ERROR;
  }

  if (parse.nonOptionsCount() != 1) {
    option::printUsage(cout, usage);
    return ERROR;
  }

  bool unknown_options = false;
  for (option::Option * opt = options[UNKNOWN_OPTION]; opt; opt = opt->next()) {
    unknown_options = true;
  }

  if (unknown_options) {
    option::printUsage(cout, usage);
    return ERROR;
  }

  // try-catch block used to detect incompatible options.
  try {
    for (int i = 0; i < parse.optionsCount(); ++i) {
      option::Option & opt = buffer[i];
      switch (opt.index()) {
        case HELP:
          // not possible, because handled further above and exits the program
        case ENGINE: engine_ = to_engine(opt.arg); break;
        case BOUND: bound_ = atoi(opt.arg); break;
        case PROP: prop_idx_ = atoi(opt.arg); break;
        case VERBOSITY: verbosity_ = atoi(opt.arg); break;
        case VCDNAME:
          vcd_name_ = opt.arg;
          if (no_witness_)
            throw PonoException(
                "Options '--vcd' and '--no-witness' are incompatible.");
          break;
        case NOWITNESS:
          no_witness_ = true;
          if (!vcd_name_.empty())
            throw PonoException(
                "Options '--vcd' and '--no-witness' are incompatible.");
          break;
        case CEGPROPHARR: ceg_prophecy_arrays_ = true; break;
        case NO_CEGP_AXIOM_RED: cegp_axiom_red_ = false; break;
        case RESET: reset_name_ = opt.arg; break;
        case RESET_BND: reset_bnd_ = atoi(opt.arg); break;
        case CLK: clock_name_ = opt.arg; break;
        case UNKNOWN_OPTION:
          // not possible because Arg::Unknown returns ARG_ILLEGAL
          // which aborts the parse with an error
          break;
      }
    }
  }
  catch (PonoException & ce) {
    cout << ce.what() << endl;
    return ERROR;
  }

  filename_ = parse.nonOption(0);

  return UNKNOWN;
}

}  // namespace pono
