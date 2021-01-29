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
  RANDOM_SEED,
  VCDNAME,
  WITNESS,
  CEGPROPHARR,
  NO_CEGP_AXIOM_RED,
  STATICCOI,
  CHECK_INVAR,
  RESET,
  RESET_BND,
  CLK,
  SMT_SOLVER,
  NO_IC3_PREGEN,
  NO_IC3_INDGEN,
  IC3_RESET_INTERVAL,
  IC3_GEN_MAX_ITER,
  IC3_FUNCTIONAL_PREIMAGE,
  NO_IC3_UNSATCORE_GEN,
  MBIC3_INDGEN_MODE,
  PROFILING_LOG_FILENAME,
  PSEUDO_INIT_PROP,
  ASSUME_PROP,
  CEGP_ABS_VALS,
  CEGP_ABS_VALS_CUTOFF,
  PROMOTE_INPUTVARS,
  SYGUS_OP_LVL
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
    "interp, mbic3, ic3ia, msat-ic3ia, sygus-pdr]." },
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
  { RANDOM_SEED,
    0,
    "",
    "random-seed",
    Arg::Numeric,
    "  --random-seed, \tRandom seed." },
  { VCDNAME,
    0,
    "",
    "vcd",
    Arg::NonEmpty,
    "  --vcd \tName of Value Change Dump (VCD) if witness exists." },
  { SMT_SOLVER,
    0,
    "",
    "smt-solver",
    Arg::NonEmpty,
    "  --smt-solver \tSMT Solver to use: btor or msat or cvc4." },
  { WITNESS,
    0,
    "",
    "witness",
    Arg::None,
    "  --witness \tPrint witness if the property is false." },
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
  { STATICCOI,
    0,
    "",
    "static-coi",
    Arg::None,
    "  --static-coi \tApply static (i.e., one-time before solving) "
    "cone-of-influence analysis." },
  { CHECK_INVAR,
    0,
    "",
    "check-invar",
    Arg::None,
    "  --check-invar \tFor engines that produce invariants, check that they "
    "hold." },
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
  { NO_IC3_PREGEN,
    0,
    "",
    "ic3-no-pregen",
    Arg::None,
    "  --ic3-no-pregen \tDisable preimage generalization in ic3." },
  { NO_IC3_INDGEN,
    0,
    "",
    "ic3-no-indgen",
    Arg::None,
    "  --ic3-no-indgen \tDisable inductive generalization in ic3." },
  { IC3_RESET_INTERVAL,
    0,
    "",
    "ic3-reset-interval",
    Arg::Numeric,
    "  --ic3-reset-interval \tNumber of check-sat queries before "
    "resetting the solver. "
    "Setting it to 0 means an unbounded number of iterations."
    "Note: some solvers don't support resetting assertions, in which "
    "case it will just fail to reset and not try again. This will be "
    "printed at verbosity 1." },
  { IC3_GEN_MAX_ITER,
    0,
    "",
    "ic3-gen-max-iter",
    Arg::Numeric,
    "  --ic3-gen-max-iter \tMax number of iterations "
    "(greater than zero) for unsatcore-based ic3 generalization. "
    "Setting it to 0 means an unbounded number of iterations." },
  { IC3_FUNCTIONAL_PREIMAGE,
    0,
    "",
    "ic3-functional-preimage",
    Arg::None,
    "  --ic3-functional-preimage \tUse functional preimage in ic3." },
  { NO_IC3_UNSATCORE_GEN,
    0,
    "",
    "no-ic3-unsatcore-gen",
    Arg::None,
    "  --no-ic3-unsatcore-gen \tDisable unsat core generalization during"
    " relative induction check. That extra generalization helps several IC3"
    " variants but also runs the risk of myopic over-generalization. Some IC3"
    " variants have better inductive generalization and do better with this"
    " option." },
  { MBIC3_INDGEN_MODE,
    0,
    "",
    "mbic3-indgen-mode",
    Arg::Numeric,
    "  --mbic3-indgen-mode \tModelBasedIC3 inductive generalization mode "
    "[0,2].\n\t"
    "0 - normal, 1 - embedded init constraint, 2 - interpolation." },
  { PROFILING_LOG_FILENAME,
    0,
    "",
    "profiling-log",
    Arg::NonEmpty,
    "  --profiling-log \tName of logfile for profiling output"
    " (requires build with linked profiling library 'gperftools')." },
  { PSEUDO_INIT_PROP,
    0,
    "",
    "pseudo-init-prop",
    Arg::None,
    "  --pseudo-init-prop \tReplace init and prop with state variables -- can "
    "extend trace by up to two steps. Recommended for use with ic3ia. "
    "Important note: will promote system to be relational" },
  { ASSUME_PROP,
    0,
    "",
    "assume-prop",
    Arg::None,
    "  --assume-prop \tenable assuming property in pre-state (default "
    "disabled)" },
  { CEGP_ABS_VALS,
    0,
    "",
    "cegp-abs-vals",
    Arg::None,
    "  --cegp-abs-vals \tabstract values in ceg-prophecy-arrays (only "
    "supported for IC3IA)" },
  { CEGP_ABS_VALS_CUTOFF,
    0,
    "",
    "cegp-abs-vals-cutoff",
    Arg::Numeric,
    "  --cegp-abs-vals-cutoff \tcutoff value for what to abstract - must be "
    "positive (default: 100)" },
  { PROMOTE_INPUTVARS,
    0,
    "",
    "promote-inputvars",
    Arg::None,
    "  --promote-inputvars \tpromote all input variables to state variables" },
  { SYGUS_OP_LVL,
      0,
      "",
      "op-lv",
      Arg::Numeric,
      "  --op-lv \toperator abstraction level (0-2, default:0) (only "
      "supported for SYGUS PDR)" },
  { 0, 0, 0, 0, 0, 0 }
};
/*********************************** end Option Handling setup
 * ***************************************/

namespace pono {

const std::string PonoOptions::default_profiling_log_filename_ = "";

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
        case RANDOM_SEED: random_seed_ = atoi(opt.arg); break;
        case VCDNAME:
          vcd_name_ = opt.arg;
          witness_ = true;  // implicitly enabling witness
          break;
        case SMT_SOLVER: {
          if (opt.arg == std::string("btor")) {
            smt_solver_ = smt::BTOR;
          } else if (opt.arg == std::string("cvc4")) {
            smt_solver_ = smt::CVC4;
          } else if (opt.arg == std::string("msat")) {
            smt_solver_ = smt::MSAT;
          } else {
            throw PonoException("Unknown solver: " + std::string(opt.arg));
            break;
          }
          break;
        }
        case WITNESS: witness_ = true; break;
        case CEGPROPHARR: ceg_prophecy_arrays_ = true; break;
        case NO_CEGP_AXIOM_RED: cegp_axiom_red_ = false; break;
        case STATICCOI: static_coi_ = true; break;
        case CHECK_INVAR: check_invar_ = true; break;
        case RESET: reset_name_ = opt.arg; break;
        case RESET_BND: reset_bnd_ = atoi(opt.arg); break;
        case CLK: clock_name_ = opt.arg; break;
        case NO_IC3_PREGEN: ic3_pregen_ = false; break;
        case NO_IC3_INDGEN: ic3_indgen_ = false; break;
        case IC3_RESET_INTERVAL: ic3_reset_interval_ = atoi(opt.arg); break;
        case IC3_GEN_MAX_ITER: ic3_gen_max_iter_ = atoi(opt.arg); break;
        case MBIC3_INDGEN_MODE:
          mbic3_indgen_mode = atoi(opt.arg);
          if (!(mbic3_indgen_mode >= 0 && mbic3_indgen_mode <= 2))
            throw PonoException(
                "--ic3-indgen-mode value must be between 0 and 2.");
          break;
        case IC3_FUNCTIONAL_PREIMAGE: ic3_functional_preimage_ = true; break;
        case NO_IC3_UNSATCORE_GEN: ic3_unsatcore_gen_ = false; break;
        case PROFILING_LOG_FILENAME:
#ifndef WITH_PROFILING
          throw PonoException(
              "Profiling requires linking to gperftools library. "
              "Please reconfigure Pono with './configure --with-profiling'.");
#else
          profiling_log_filename_ = opt.arg;
#endif
          break;
        case PSEUDO_INIT_PROP: pseudo_init_prop_ = true; break;
        case ASSUME_PROP: assume_prop_ = true; break;
        case CEGP_ABS_VALS: cegp_abs_vals_ = true; break;
        case CEGP_ABS_VALS_CUTOFF: cegp_abs_vals_cutoff_ = atoi(opt.arg); break;
        case PROMOTE_INPUTVARS: promote_inputvars_ = true; break;
        case SYGUS_OP_LVL: sygus_use_operator_abstraction_ = atoi(opt.arg); break;
        case UNKNOWN_OPTION:
          // not possible because Arg::Unknown returns ARG_ILLEGAL
          // which aborts the parse with an error
          break;
      }
    }

    if (smt_solver_ != smt::MSAT && engine_ == Engine::INTERP) {
      throw PonoException(
          "Interpolation engine can be only used with '--smt-solver msat'.");
    }

    if (ceg_prophecy_arrays_ && smt_solver_ != smt::MSAT) {
      throw PonoException(
          "Counterexample-guided prophecy only supported with MathSAT so far");
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
