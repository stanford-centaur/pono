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

#include <climits>
#include <iostream>
#include <string>
#include <vector>

#include "optionparser.h"
#include "utils/exceptions.h"

#ifndef PONO_VERSION
#define PONO_VERSION "unknown"
#endif

using namespace std;

/************************************* Option Handling setup
 * *****************************************/
// from optionparser-1.7 examples -- example_arg.cc
enum optionIndex
{
  UNKNOWN_OPTION,
  HELP,
  VERSION,
  ENGINE,
  BOUND,
  PROP,
  VERBOSITY,
  RANDOM_SEED,
  VCDNAME,
  WITNESS,
  BTOR2_WITNESS_NAME,
  JUSTICE,
  JUSTICE_TRANSLATOR,
  STATICCOI,
  SHOW_INVAR,
  CHECK_INVAR,
  RESET,
  RESET_BND,
  CLK,
  SMT_SOLVER,
  SMT_INTERPOLATOR,
  LOGGING_SMT_SOLVER,
  PRINTING_SMT_SOLVER,
  NO_IC3_PREGEN,
  NO_IC3_INDGEN,
  IC3_GEN_MAX_ITER,
  IC3_FUNCTIONAL_PREIMAGE,
  NO_IC3_UNSATCORE_GEN,
  NO_IC3IA_REDUCE_PREDS,
  NO_IC3IA_TRACK_IMPORTANT_VARS,
  NO_IC3IA_SIM_CEX,
  NO_IC3SA_FUNC_REFINE,
  MBIC3_INDGEN_MODE,
  PROFILING_LOG_FILENAME,
  PSEUDO_INIT_PROP,
  ASSUME_PROP,
  CEGPROPHARR,
  NO_CEGP_TIMED_AXIOM_RED,
  NO_CEGP_CONSEC_AXIOM_RED,
  NO_CEGP_NONCONSEC_AXIOM_RED,
  CEGP_FORCE_RESTART,
  CEGP_ABS_VALS,
  CEGP_ABS_VALS_CUTOFF,
  CEGP_STRONG_ABSTRACTION,
  CEG_BV_ARITH,
  CEG_BV_ARITH_MIN_BW,
  PROMOTE_INPUTVARS,
  SYGUS_OP_LVL,
  SYGUS_TERM_MODE,
  IC3SA_INITIAL_TERMS_LVL,
  IC3SA_INTERP,
  PRINT_WALL_TIME,
  BMC_BOUND_START,
  BMC_BOUND_STEP,
  BMC_NEG_INIT_STEP,
  BMC_EXPONENTIAL_STEP,
  BMC_SINGLE_BAD_STATE,
  BMC_NEG_BAD_STEP,
  BMC_MIN_CEX_LIN_SEARCH,
  BMC_MIN_CEX_LESS_INC_BIN_SEARCH,
  BMC_NEG_BAD_STEP_ALL,
  BMC_ALLOW_NON_MINIMAL_CEX,
  KIND_NO_SIMPLE_PATH_CHECK,
  KIND_EAGER_SIMPLE_PATH_CHECK,
  KIND_NO_MULTI_CALL_SIMPLE_PATH_CHECK,
  KIND_NO_IND_CHECK_INIT_STATES,
  KIND_NO_IND_CHECK,
  KIND_NO_IND_CHECK_PROPERTY,
  KIND_ONE_TIME_BASE_CHECK,
  KIND_BOUND_STEP,
  NO_INTERP_FRONTIER_SIMPL,
  INTERP_PROPS,
  INTERP_EAGER_UNROLL,
  INTERP_BACKWARD
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
  { VERSION,
    0,
    "",
    "version",
    Arg::None,
    "  --version \tPrint version and exit." },
  { ENGINE,
    0,
    "e",
    "engine",
    Arg::NonEmpty,
    "  --engine, -e <engine> \tSelect engine from [bmc, bmc-sp, ind, interp, "
    "ismc, mbic3, ic3bits, ic3ia, msat-ic3ia, ic3sa, sygus-pdr]." },
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
    "  --verbosity, -v \tVerbosity for printing to stderr." },
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
    "  --vcd <filename> \tName of Value Change Dump (VCD) if witness exists." },
  { SMT_SOLVER,
    0,
    "",
    "smt-solver",
    Arg::NonEmpty,
    "  --smt-solver \tSMT Solver to use: btor, bzla, msat, yices2, or cvc5. "
    "(default: bzla)" },
  { SMT_INTERPOLATOR,
    0,
    "",
    "smt-interpolator",
    Arg::NonEmpty,
    "  --smt-interpolator \tSMT Solver used for interpolation: msat or cvc5. "
#ifdef WITH_MSAT
    "(default: msat)"
#else
    "(default: cvc5)"
#endif
  },
  { LOGGING_SMT_SOLVER,
    0,
    "",
    "logging-smt-solver",
    Arg::None,
    "  --logging-smt-solver \tUse Smt-Switch logging solver which guarantees "
    "the exact term structure that was created. Good for avoiding term "
    "rewriting at the API level or sort aliasing. (default: false)" },
  { PRINTING_SMT_SOLVER,
    0,
    "",
    "printing-smt-solver",
    Arg::None,
    "  --printing-smt-solver \tDump all SMT queries to standard error output "
    "in SMT-LIB format while solving. Uses smt-switch's create_printing_solver "
    "function. (default: false)" },
  { WITNESS,
    0,
    "",
    "witness",
    Arg::None,
    "  --witness \tPrint witness if the property is false." },
  { BTOR2_WITNESS_NAME,
    0,
    "",
    "dump-btor2-witness",
    Arg::NonEmpty,
    "  --dump-btor2-witness <filename> \t"
    "Dump the Btor2 witness into the specified file." },
  { JUSTICE,
    0,
    "j",
    "justice",
    Arg::None,
    "  --justice, -j \tCheck justice property, instead of safety, at the given "
    "index. (default: false)" },
  { JUSTICE_TRANSLATOR,
    0,
    "",
    "justice-translator",
    Arg::NonEmpty,
    "  --justice-translator <algorithm> \tSelect liveness to safety "
    "translation algorithm from [l2s]." },
  { STATICCOI,
    0,
    "",
    "static-coi",
    Arg::None,
    "  --static-coi \tApply static (i.e., one-time before solving) "
    "cone-of-influence analysis." },
  { SHOW_INVAR,
    0,
    "",
    "show-invar",
    Arg::None,
    "  --show-invar \tFor engines that produce invariants, show the "
    "invariant" },
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
    "~ for negative reset)" },
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
    "supports starting at 0 and toggling each step)" },
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
  { IC3_GEN_MAX_ITER,
    0,
    "",
    "ic3-gen-max-iter",
    Arg::Numeric,
    "  --ic3-gen-max-iter \tMax number of iterations (greater than zero) for "
    "unsatcore-based ic3 generalization. Setting it to 0 means an unbounded "
    "number of iterations." },
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
    "  --no-ic3-unsatcore-gen \tDisable unsat core generalization during "
    "relative induction check. That extra generalization helps several IC3 "
    "variants but also runs the risk of myopic over-generalization. Some IC3 "
    "variants have better inductive generalization and do better with this "
    "option." },
  { NO_IC3IA_REDUCE_PREDS,
    0,
    "",
    "no-ic3ia-reduce-preds",
    Arg::None,
    "  --no-ic3ia-reduce-preds \tDisable unsat core based predicate "
    "minimization" },
  { NO_IC3IA_TRACK_IMPORTANT_VARS,
    0,
    "",
    "no-ic3ia-track-important-vars",
    Arg::None,
    "  --no-ic3ia-track-important-vars \tIgnore tracked important variables "
    "when picking predicates." },
  { NO_IC3IA_SIM_CEX,
    0,
    "",
    "no-ic3ia-sim-cex",
    Arg::None,
    "  --no-ic3ia-sim-cex \tDo not simulate abstract cex during refinement, "
    "perform BMC instead." },
  { NO_IC3SA_FUNC_REFINE,
    0,
    "",
    "no-ic3sa-func-refine",
    Arg::None,
    "  --no-ic3sa-func-refine \tDisable functional unrolling attempt in "
    "IC3SA." },
  { MBIC3_INDGEN_MODE,
    0,
    "",
    "mbic3-indgen-mode",
    Arg::Numeric,
    "  --mbic3-indgen-mode \tModelBasedIC3 inductive generalization mode "
    "[0,2].\n\t0 - normal, 1 - embedded init constraint, 2 - interpolation." },
  { PROFILING_LOG_FILENAME,
    0,
    "",
    "profiling-log",
    Arg::NonEmpty,
    "  --profiling-log \tName of logfile for profiling output (requires build "
    "with linked profiling library 'gperftools')." },
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
  { CEGPROPHARR,
    0,
    "",
    "ceg-prophecy-arrays",
    Arg::None,
    "  --ceg-prophecy-arrays \tUse counter-example guided prophecy for "
    "arrays." },
  { NO_CEGP_TIMED_AXIOM_RED,
    0,
    "",
    "no-cegp-timed-axiom-red",
    Arg::None,
    "  --no-cegp-timed-axiom-red \tDon't reduce enumerated axioms in "
    "CEG-Prophecy with unsat cores." },
  { NO_CEGP_CONSEC_AXIOM_RED,
    0,
    "",
    "no-cegp-consec-axiom-red",
    Arg::None,
    "  --no-cegp-consec-axiom-red \tDon't reduce consecutive axioms in "
    "CEG-Prophecy before adding to transition system." },
  { NO_CEGP_NONCONSEC_AXIOM_RED,
    0,
    "",
    "no-cegp-nonconsec-axiom-red",
    Arg::None,
    "  --no-cegp-nonconsec-axiom-red \tDon't reduce non-consecutive axioms in "
    "CEG-Prophecy before creating prophecy variables." },
  { CEGP_FORCE_RESTART,
    0,
    "",
    "cegp-force-restart",
    Arg::None,
    "  --cegp-force-restart \tForce underlying engine to restart after "
    "refinement." },
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
  { CEGP_STRONG_ABSTRACTION,
    0,
    "",
    "cegp-strong-abs",
    Arg::None,
    "  --cegp-strong-abs \tUse strong abstraction in CEGP -- no equality UFs "
    "(default: false)" },
  { CEG_BV_ARITH,
    0,
    "",
    "ceg-bv-arith",
    Arg::None,
    "  --ceg-bv-arith \tabstraction-refinement for the BV arithmetic operators "
    "(mul, div, rem, mod)" },
  { CEG_BV_ARITH_MIN_BW,
    0,
    "",
    "ceg-bv-arith-min-bw",
    Arg::Numeric,
    "  --ceg-bv-arith-min-bw \tminimum bitwidth of operators to abstract - "
    "must be positive (default: 16) " },
  { PROMOTE_INPUTVARS,
    0,
    "",
    "promote-inputvars",
    Arg::None,
    "  --promote-inputvars \tpromote all input variables to state variables" },
  { SYGUS_OP_LVL,
    0,
    "",
    "sygus-op-lv",
    Arg::Numeric,
    "  --sygus-op-lv \toperator abstraction level (0-2, default:0) (only for "
    "SYGUS PDR)" },
  { SYGUS_TERM_MODE,
    0,
    "",
    "sygus-term-mode",
    Arg::Numeric,
    "  --sygus-term-mode \tterm generation mode (0-4, default: 4 AUTO) (0: "
    "more replace, 1: v/c ext 2: v/c split 3: v/c lt/le )" },
  { IC3SA_INITIAL_TERMS_LVL,
    0,
    "",
    "ic3sa-initial-terms-lvl",
    Arg::Numeric,
    "  --ic3sa-initial-terms-lvl \tConfigures where to find terms for the "
    "initial abstraction. Higher numbers means more terms and predicates will "
    "be included in the initial abstraction [0-4] (default: 4)." },
  { IC3SA_INTERP,
    0,
    "",
    "ic3sa-interp",
    Arg::None,
    "  --ic3sa-interp \tuse interpolants to find more terms during refinement "
    "(default: off)" },
  { PRINT_WALL_TIME,
    0,
    "",
    "print-wall-time",
    Arg::None,
    "  --print-wall-time \tPrint wall clock time of entire execution" },
  { BMC_BOUND_START,
    0,
    "",
    "bmc-bound-start",
    Arg::Numeric,
    "  --bmc-bound-start \tBound (unrolling depth) to start cex search in BMC "
    "(default: 0)" },
  { BMC_BOUND_STEP,
    0,
    "",
    "bmc-bound-step",
    Arg::Numeric,
    "  --bmc-bound-step \tAmount by which bound (unrolling depth) for cex "
    "search is increased in BMC (default: 1). For values greater than 1, BMC "
    "searches for cex in intervals of size '--bmc-bound-step'." },
  { BMC_NEG_INIT_STEP,
    0,
    "",
    "bmc-neg-init-step",
    Arg::None,
    "  --bmc-neg-init-step \tAdd negated initial state constraint in BMC steps "
    "k > 0 (default: false)." },
  { BMC_EXPONENTIAL_STEP,
    0,
    "",
    "bmc-exponential-step",
    Arg::None,
    "  --bmc-exponential-step \tDouble BMC bound in each step starting at "
    "'bmc-bound-start' (default: false, explores bounds 0, 1, 2, 4,...)." },
  { BMC_SINGLE_BAD_STATE,
    0,
    "",
    "bmc-single-bad-state",
    Arg::None,
    "  --bmc-single-bad-state \tEXPERT OPTION: add a single bad state literal "
    "for current bound k rather than a disjunctive term covering the checked "
    "interval; counterexamples may be missed. (default: false)." },
  { BMC_NEG_BAD_STEP,
    0,
    "",
    "bmc-neg-bad-step",
    Arg::None,
    "  --bmc-neg-bad-step \tAdd negated bad state constraint in BMC steps k > "
    "0 (default: false)." },
  { BMC_NEG_BAD_STEP_ALL,
    0,
    "",
    "bmc-neg-bad-step-all",
    Arg::None,
    "  --bmc-neg-bad-step-all \tEXPERT OPTION: like '--bmc-neg-bad-step' but "
    "add negated bad state constraint in ALL BMC steps k > 0 (default: false). "
    "When combined with --bmc-single-bad-state, this option may cause "
    "overconstraining the problem in certain corner cases." },
  { BMC_MIN_CEX_LIN_SEARCH,
    0,
    "",
    "bmc-min-cex-linear-search",
    Arg::None,
    "  --bmc-min-cex-linear-search \tApply linear instead of binary search for "
    "minimal cex after a cex was found in current interval" },
  { BMC_MIN_CEX_LESS_INC_BIN_SEARCH,
    0,
    "",
    "bmc-min-cex-less-inc-bin-search",
    Arg::None,
    "  --bmc-min-cex-less-inc-bin-search \tApply less incremental variant of "
    "binary search for minimal cex after a cex was found in current interval" },
  { BMC_ALLOW_NON_MINIMAL_CEX,
    0,
    "",
    "bmc-allow-non-minimal-cex",
    Arg::None,
    "  --bmc-allow-non-minimal-cex \tDo not search for minimal cex within an "
    "interval; instead, terminate immediately (reported bound of cex is an "
    "upper bound of actual cex)" },
  { KIND_NO_SIMPLE_PATH_CHECK,
    0,
    "",
    "kind-no-simple-path-check",
    Arg::None,
    "  --kind-no-simple-path-check \tSkip simple path check in k-induction "
    "(WARNING: might cause incompleteness)" },
  { KIND_EAGER_SIMPLE_PATH_CHECK,
    0,
    "",
    "kind-eager-simple-path-check",
    Arg::None,
    "  --kind-eager-simple-path-check \tEager simple path check in k-induction "
    "(default: lazy check)" },
  { KIND_NO_MULTI_CALL_SIMPLE_PATH_CHECK,
    0,
    "",
    "kind-no-multi-call-simple-path-check",
    Arg::None,
    "  --kind-no-multi-call-simple-path-check \tTry to avoid multiple solver "
    "calls in lazy simple path check in k-induction" },
  { KIND_NO_IND_CHECK_INIT_STATES,
    0,
    "",
    "kind-no-ind-check-init-states",
    Arg::None,
    "  --kind-no-ind-check-init-states \tK-induction: skip checking inductive "
    "case based on initial states" },
  { KIND_NO_IND_CHECK,
    0,
    "",
    "kind-no-ind-check",
    Arg::None,
    "  --kind-no-ind-check \tK-induction: skip inductive case checks; implies "
    "'--kind-no-ind-check-init-states' and '--kind-no-ind-check-property' "
    "(WARNING: will cause incompleteness on most problem instances)" },
  { KIND_NO_IND_CHECK_PROPERTY,
    0,
    "",
    "kind-no-ind-check-property",
    Arg::None,
    "  --kind-no-ind-check-property \tK-induction: skip checking inductive "
    "case based on property (WARNING: will cause incompleteness on most "
    "problem instances)" },
  { KIND_ONE_TIME_BASE_CHECK,
    0,
    "",
    "kind-one-time-base-check",
    Arg::None,
    "  --kind-one-time-base-check \tK-induction: check base case only once "
    "after inductive check was unsatisfiable (WARNING: counterexamples might "
    "be missed)" },
  { KIND_BOUND_STEP,
    0,
    "",
    "kind-bound-step",
    Arg::Numeric,
    "  --kind-bound-step \tAmount by which bound (unrolling depth) is "
    "increased in k-induction (default: 1)" },
  { NO_INTERP_FRONTIER_SIMPL,
    0,
    "",
    "no-interp-frontier-simpl",
    Arg::None,
    "  --no-interp-frontier-simpl \tDisable frontier-set simplification in "
    "interp engine" },
  { INTERP_PROPS,
    0,
    "",
    "interp-props",
    Arg::NonEmpty,
    "  --interp-props \tSpecifies at which time frames properties are "
    "considered when computing interpolants: all (default) and first-and-last "
    "(WARNING: choosing 'fist-and-last' could cause incompleteness on some "
    "instances)" },
  { INTERP_EAGER_UNROLL,
    0,
    "",
    "interp-eager-unroll",
    Arg::None,
    "  --interp-eager-unroll \tUnroll the transition system eagerly in interp "
    "engine" },
  { INTERP_BACKWARD,
    0,
    "",
    "interp-backward",
    Arg::None,
    "  --interp-backward \tCompute interpolants in a backward manner, "
    "i.e., not(itp(B, A)), in interp engine "
    "(forward, i.e., itp(A, B), if not specified)" },
  { 0, 0, 0, 0, 0, 0 }
};
/*********************************** end Option Handling setup
 * ***************************************/

namespace pono {

const std::unordered_set<Engine> ic3_variants_set({ IC3_BOOL,
                                                    IC3_BITS,
                                                    MBIC3,
                                                    IC3IA_ENGINE,
                                                    MSAT_IC3IA,
                                                    IC3SA_ENGINE,
                                                    SYGUS_PDR });

const std::unordered_set<Engine> & ic3_variants() { return ic3_variants_set; }

const std::string PonoOptions::default_profiling_log_filename_ = "";

Engine PonoOptions::to_engine(std::string s)
{
  if (str2engine.find(s) != str2engine.end()) {
    return str2engine.at(s);
  } else {
    throw PonoException("Unrecognized engine: " + s);
  }
}

JusticeTranslator PonoOptions::to_justice_translator(std::string s)
{
  if (str2livenessalg.find(s) != str2livenessalg.end()) {
    return str2livenessalg.at(s);
  } else {
    throw PonoException("Unrecognized algorithm: " + s);
  }
}

// Parse command line options given by 'argc' and 'argv' and set
// respective options in the 'pono_options' object.
// Returns 'ERROR' if there is something wrong with the given options
// or 'UNKNOWN' instead.
ProverResult PonoOptions::parse_and_set_options(int argc,
                                                char ** argv,
                                                bool expect_file)
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

  if (options[VERSION]) {
    cout << PONO_VERSION << endl;
    return ERROR;
  }

  if (expect_file && parse.nonOptionsCount() != 1
      || parse.nonOptionsCount() > 1) {
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
        case BOUND:
          bound_ = std::stoul(opt.arg);
          if (bound_ >= INT_MAX) {
            throw PonoException("--bound must be less than "
                                + std::to_string(INT_MAX) + ".");
          }
          break;
        case PROP: prop_idx_ = std::stoul(opt.arg); break;
        case VERBOSITY: verbosity_ = std::stoul(opt.arg); break;
        case RANDOM_SEED: random_seed_ = std::stoul(opt.arg); break;
        case VCDNAME:
          vcd_name_ = opt.arg;
          witness_ = true;  // implicitly enabling witness
          break;
        case BTOR2_WITNESS_NAME:
          btor2_witness_name_ = opt.arg;
          witness_ = true;  // implicitly enabling witness
          break;
        case SMT_SOLVER: {
          if (opt.arg == std::string("btor")) {
            smt_solver_ = smt::BTOR;
          } else if (opt.arg == std::string("bzla")) {
            smt_solver_ = smt::BZLA;
          } else if (opt.arg == std::string("cvc5")) {
            smt_solver_ = smt::CVC5;
          } else if (opt.arg == std::string("msat")) {
            smt_solver_ = smt::MSAT;
          } else if (opt.arg == std::string("yices2")) {
            smt_solver_ = smt::YICES2;
          } else {
            throw PonoException("Unknown solver: " + std::string(opt.arg));
            break;
          }
          break;
        }
        case SMT_INTERPOLATOR: {
          if (opt.arg == std::string("cvc5")) {
            smt_interpolator_ = smt::CVC5_INTERPOLATOR;
          } else if (opt.arg == std::string("msat")) {
            smt_interpolator_ = smt::MSAT_INTERPOLATOR;
          } else {
            throw PonoException("Unknown interpolator: "
                                + std::string(opt.arg));
            break;
          }
          break;
        }
        case LOGGING_SMT_SOLVER: logging_smt_solver_ = true; break;
        case PRINTING_SMT_SOLVER: printing_smt_solver_ = true; break;
        case WITNESS: witness_ = true; break;
        case JUSTICE: justice_ = true; break;
        case JUSTICE_TRANSLATOR:
          justice_translator_ = to_justice_translator(opt.arg);
          break;
        case STATICCOI: static_coi_ = true; break;
        case SHOW_INVAR: show_invar_ = true; break;
        case CHECK_INVAR: check_invar_ = true; break;
        case RESET: reset_name_ = opt.arg; break;
        case RESET_BND: reset_bnd_ = std::stoul(opt.arg); break;
        case CLK: clock_name_ = opt.arg; break;
        case NO_IC3_PREGEN: ic3_pregen_ = false; break;
        case NO_IC3_INDGEN: ic3_indgen_ = false; break;
        case IC3_GEN_MAX_ITER: ic3_gen_max_iter_ = std::stoul(opt.arg); break;
        case MBIC3_INDGEN_MODE:
          mbic3_indgen_mode = std::stoul(opt.arg);
          if (!(mbic3_indgen_mode >= 0 && mbic3_indgen_mode <= 2))
            throw PonoException(
                "--ic3-indgen-mode value must be between 0 and 2.");
          break;
        case IC3_FUNCTIONAL_PREIMAGE: ic3_functional_preimage_ = true; break;
        case NO_IC3_UNSATCORE_GEN: ic3_unsatcore_gen_ = false; break;
        case NO_IC3IA_REDUCE_PREDS: ic3ia_reduce_preds_ = false;
        case NO_IC3IA_TRACK_IMPORTANT_VARS: ic3ia_track_important_vars_ = false;
        case NO_IC3IA_SIM_CEX: ic3ia_sim_cex_ = false; break;
        case NO_IC3SA_FUNC_REFINE: ic3sa_func_refine_ = false; break;
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
        case CEGPROPHARR: ceg_prophecy_arrays_ = true; break;
        case NO_CEGP_TIMED_AXIOM_RED: cegp_timed_axiom_red_ = false; break;
        case NO_CEGP_CONSEC_AXIOM_RED: cegp_consec_axiom_red_ = false; break;
        case NO_CEGP_NONCONSEC_AXIOM_RED:
          cegp_nonconsec_axiom_red_ = false;
          break;
        case CEGP_FORCE_RESTART: cegp_force_restart_ = true; break;
        case CEGP_ABS_VALS: cegp_abs_vals_ = true; break;
        case CEGP_ABS_VALS_CUTOFF:
          cegp_abs_vals_cutoff_ = std::stoul(opt.arg);
          break;
        case CEGP_STRONG_ABSTRACTION: cegp_strong_abstraction_ = true; break;
        case CEG_BV_ARITH: ceg_bv_arith_ = true; break;
        case CEG_BV_ARITH_MIN_BW:
          ceg_bv_arith_min_bw_ = std::stoul(opt.arg);
          break;
        case PROMOTE_INPUTVARS: promote_inputvars_ = true; break;
        case SYGUS_OP_LVL:
          sygus_use_operator_abstraction_ = std::stoul(opt.arg);
          break;
        case SYGUS_TERM_MODE:
          sygus_term_mode_ = SyGuSTermMode(std::stoul(opt.arg));
          break;
        case IC3SA_INITIAL_TERMS_LVL: {
          ic3sa_initial_terms_lvl_ = std::stoul(opt.arg);
          if (ic3sa_initial_terms_lvl_ > 4) {
            throw PonoException(
                "--ic3sa-initial-terms-lvl must be an integer in [0, 4]");
          }
          break;
        }
        case IC3SA_INTERP: ic3sa_interp_ = true; break;
        case PRINT_WALL_TIME: print_wall_time_ = true; break;
        case BMC_BOUND_START: bmc_bound_start_ = std::stoul(opt.arg); break;
        case BMC_BOUND_STEP:
          bmc_bound_step_ = std::stoul(opt.arg);
          if (bmc_bound_step_ == 0)
            throw PonoException("--bmc-bound-step must be greater than 0");
          break;
        case BMC_NEG_INIT_STEP: bmc_neg_init_step_ = true; break;
        case BMC_EXPONENTIAL_STEP: bmc_exponential_step_ = true; break;
        case BMC_SINGLE_BAD_STATE: bmc_single_bad_state_ = true; break;
        case BMC_NEG_BAD_STEP:
          bmc_neg_bad_step_ = true;
          if (bmc_neg_bad_step_all_)
            throw PonoException(
                "--bmc-neg-bad-step-all cannot be combined with "
                "'--bmc-neg-bad-step'");
          break;
        case BMC_NEG_BAD_STEP_ALL:
          bmc_neg_bad_step_all_ = true;
          if (bmc_neg_bad_step_)
            throw PonoException(
                "--bmc-neg-bad-step cannot be combined with "
                "'--bmc-neg-bad-step-all'");
          break;
        case BMC_MIN_CEX_LIN_SEARCH: bmc_min_cex_linear_search_ = true; break;
        case BMC_MIN_CEX_LESS_INC_BIN_SEARCH:
          bmc_min_cex_less_inc_bin_search_ = true;
          break;
        case BMC_ALLOW_NON_MINIMAL_CEX:
          bmc_allow_non_minimal_cex_ = true;
          break;
        case KIND_NO_SIMPLE_PATH_CHECK:
          kind_no_simple_path_check_ = true;
          break;
        case KIND_EAGER_SIMPLE_PATH_CHECK:
          kind_eager_simple_path_check_ = true;
          break;
        case KIND_NO_MULTI_CALL_SIMPLE_PATH_CHECK:
          kind_no_multi_call_simple_path_check_ = true;
          break;
        case KIND_NO_IND_CHECK_INIT_STATES:
          kind_no_ind_check_init_states_ = true;
          break;
        case KIND_NO_IND_CHECK:
          kind_no_ind_check_ = true;
          kind_no_ind_check_init_states_ = true;
          kind_no_ind_check_property_ = true;
          break;
        case KIND_NO_IND_CHECK_PROPERTY:
          kind_no_ind_check_property_ = true;
          break;
        case KIND_ONE_TIME_BASE_CHECK: kind_one_time_base_check_ = true; break;
        case KIND_BOUND_STEP:
          kind_bound_step_ = std::stoul(opt.arg);
          if (kind_bound_step_ == 0)
            throw PonoException("--kind-bound-step must be greater than 0");
          break;
        case NO_INTERP_FRONTIER_SIMPL:
          interp_frontier_set_simpl_ = false;
          break;
        case INTERP_PROPS:
          if (opt.arg == std::string("all")) {
            interp_props_ = InterpPropsEnum::INTERP_ALL_PROPS;
          } else if (opt.arg == std::string("first-and-last")) {
            interp_props_ = InterpPropsEnum::INTERP_FIRST_AND_LAST_PROPS;
          } else {
            throw PonoException("Unknown --interp-props option: "
                                + std::string(opt.arg));
          }
          break;
        case INTERP_EAGER_UNROLL: interp_eager_unroll_ = true; break;
        case INTERP_BACKWARD: interp_backward_ = true; break;
        case UNKNOWN_OPTION:
          // not possible because Arg::Unknown returns ARG_ILLEGAL
          // which aborts the parse with an error
          break;
      }
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

  if (expect_file) {
    filename_ = parse.nonOption(0);
  }

  return UNKNOWN;
}

ProverResult PonoOptions::parse_and_set_options(std::vector<std::string> & opts,
                                                bool expect_file)
{
  // add one for dummy program name
  int size = opts.size() + 1;
  std::vector<char *> cstrings({ std::string("pono").data() });
  cstrings.reserve(size);
  for (auto & o : opts) {
    cstrings.push_back(o.data());
  }
  return parse_and_set_options(size, cstrings.data(), expect_file);
}

string to_string(Engine e)
{
  string res;
  switch (e) {
    case BMC: {
      res = "bmc";
      break;
    }
    case BMC_SP: {
      res = "bmc-sp";
      break;
    }
    case KIND: {
      res = "ind";
      break;
    }
    case INTERP: {
      res = "interp";
      break;
    }
    case ISMC: {
      res = "ismc";
      break;
    }
    case MBIC3: {
      res = "mbic3";
      break;
    }
    case IC3_BOOL: {
      res = "ic3bool";
      break;
    }
    case IC3_BITS: {
      res = "ic3bits";
      break;
    }
    case IC3IA_ENGINE: {
      res = "ic3ia";
      break;
    }
    case MSAT_IC3IA: {
      res = "msat-ic3ia";
      break;
    }
    case IC3SA_ENGINE: {
      res = "ic3sa";
      break;
    }
    case SYGUS_PDR: {
      res = "sygus-pdr";
      break;
    }
    default: {
      throw PonoException("Unhandled engine: " + std::to_string(e));
    }
  }
  return res;
}

string to_string(JusticeTranslator jt)
{
  string res;
  switch (jt) {
    case LIVENESS_TO_SAFETY: res = "l2s"; break;
  }
  return res;
}

ostream & operator<<(ostream & o, Engine e)
{
  o << to_string(e);
  return o;
}

ostream & operator<<(ostream & o, JusticeTranslator jt)
{
  o << to_string(jt);
  return o;
}
}  // namespace pono
