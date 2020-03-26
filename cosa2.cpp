/*********************                                                        */
/*! \file 
 ** \verbatim
 ** Top contributors (to current version):
 **   Makai Mann, Ahmed Irfan
 ** This file is part of the cosa2 project.
 ** Copyright (c) 2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief 
 **
 ** 
 **/


#include <iostream>
#include "assert.h"

#include "optionparser.h"
#include "smt-switch/boolector_factory.h"
#ifdef WITH_MSAT
  #include "smt-switch/msat_factory.h"
#endif

#include "bmc.h"
#include "bmc_simplepath.h"
#include "core/fts.h"
#include "defaults.h"
#include "frontends/btor2_encoder.h"
#include "frontends/smv_encoder.h"
#include "interpolant.h"
#include "kinduction.h"
#include "printers/btor2_witness_printer.h"
#include "prop.h"
#include "utils/logger.h"

using namespace cosa;
using namespace smt;
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
  VERBOSITY
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
    "USAGE: cosa2 [options] <btor file>\n\n"
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
  { 0, 0, 0, 0, 0, 0 }
};
/*********************************** end Option Handling setup
 * ***************************************/

ProverResult check_prop(Engine engine,
                        unsigned int bound,
                        Property & p,
                        SmtSolver & s,
                        SmtSolver & second_solver,
                        std::vector<UnorderedTermMap> & cex)
{
  logger.log(1, "Solving property: {}", p.prop());

  logger.log(3, "INIT:\n{}", p.transition_system().init());
  logger.log(3, "TRANS:\n{}", p.transition_system().trans());

  std::shared_ptr<Prover> prover;
  if (engine == BMC) {
    prover = std::make_shared<Bmc>(p, s);
  } else if (engine == BMC_SP) {
    prover = std::make_shared<BmcSimplePath>(p, s);
  } else if (engine == KIND) {
    prover = std::make_shared<KInduction>(p, s);
  } else if (engine == INTERP) {
    assert(second_solver != NULL);
    prover = std::make_shared<InterpolantMC>(p, s, second_solver);
  } else {
    throw CosaException("Unimplemented engine.");
  }

  ProverResult r = prover->check_until(bound);
  if (r == FALSE) {
    prover->witness(cex);
  }
  return r;
}

int main(int argc, char ** argv)
{
  argc -= (argc > 0);
  argv += (argc > 0);  // skip program name argv[0] if present
  option::Stats stats(usage, argc, argv);
  std::vector<option::Option> options(stats.options_max);
  std::vector<option::Option> buffer(stats.buffer_max);
  option::Parser parse(usage, argc, argv, &options[0], &buffer[0]);

  if (parse.error()) {
    return 3;
  }

  if (options[HELP] || argc == 0) {
    option::printUsage(cout, usage);
    return 2;  // unknown is 2
  }

  if (parse.nonOptionsCount() != 1) {
    option::printUsage(cout, usage);
    return 3;
  }

  bool unknown_options = false;
  for (option::Option * opt = options[UNKNOWN_OPTION]; opt; opt = opt->next()) {
    unknown_options = true;
  }

  if (unknown_options) {
    option::printUsage(cout, usage);
    return 3;
  }

  Engine engine = default_engine;
  unsigned int prop_idx = default_prop_idx;
  unsigned int bound = default_bound;
  unsigned int verbosity = default_verbosity;

  for (int i = 0; i < parse.optionsCount(); ++i) {
    option::Option & opt = buffer[i];
    switch (opt.index()) {
      case HELP:
        // not possible, because handled further above and exits the program
      case ENGINE: engine = to_engine(opt.arg); break;
      case BOUND: bound = atoi(opt.arg); break;
      case PROP: prop_idx = atoi(opt.arg); break;
      case VERBOSITY: verbosity = atoi(opt.arg); break;
      case UNKNOWN_OPTION:
        // not possible because Arg::Unknown returns ARG_ILLEGAL
        // which aborts the parse with an error
        break;
    }
  }

  // set logger verbosity -- can only be set once
  logger.set_verbosity(verbosity);

  string filename(parse.nonOption(0));

  int status_code = 3;

  try {
    SmtSolver s;
    SmtSolver second_solver;
    if (engine == INTERP) {
      #ifdef WITH_MSAT
      // need mathsat for interpolant based model checking
      s = MsatSolverFactory::create();
      second_solver = MsatSolverFactory::create_interpolating_solver();
      #else
      throw CosaException("Interpolation-based model checking requires MathSAT and "
                          "this version of cosa2 is built without MathSAT.\nPlease "
                          "setup smt-switch with MathSAT and reconfigure using --with-msat.\n"
                          "Note: MathSAT has a custom license and you must assume all "
                          "responsibility for meeting the license requirements.");
      #endif
    } else {
      // boolector is faster but doesn't support interpolants
      s = BoolectorSolverFactory::create();
      s->set_opt("produce-models", "true");
      s->set_opt("incremental", "true");
    }

    // TODO: make this less ugly, just need to keep it in scope if using
    //       it would be better to have a generic encoder
    //       and also only create the transition system once
    ProverResult r;
    string file_ext = filename.substr(filename.find_last_of(".") + 1);
    if (file_ext == "btor2" || file_ext == "btor") {
      logger.log(2, "Parsing BTOR2 file: {}", filename);

      FunctionalTransitionSystem fts(s);
      BTOR2Encoder btor_enc(filename, fts);
      const TermVec & propvec = btor_enc.propvec();
      unsigned int num_props = propvec.size();
      if (prop_idx >= num_props) {
        throw CosaException(
            "Property index " + to_string(prop_idx)
            + " is greater than the number of properties in file " + filename
            + " (" + to_string(num_props) + ")");
      }
      Term prop = propvec[prop_idx];
      Property p(fts, prop);
      vector<UnorderedTermMap> cex;
      r = check_prop(engine, bound, p, s, second_solver, cex);

      // print btor output
      if (r == FALSE) {
        cout << "sat" << endl;
        cout << "b" << prop_idx << endl;
        if (cex.size()) {
          print_witness_btor(btor_enc, cex);
        }
        status_code = 1;
      } else if (r == TRUE) {
        cout << "unsat" << endl;
        cout << "b" << prop_idx << endl;
        status_code = 0;
      } else {
        cout << "unknown" << endl;
        cout << "b" << prop_idx << endl;
        status_code = 2;
      }

    } else if (file_ext == "smv") {
      logger.log(2, "Parsing SMV file: {}", filename);
      RelationalTransitionSystem rts(s);
      smvEncoder smv_enc(filename, rts);
      const TermVec & propvec = smv_enc.propvec();
      unsigned int num_props = propvec.size();
      if (prop_idx >= num_props) {
        throw CosaException(
            "Property index " + to_string(prop_idx)
            + " is greater than the number of properties in file " + filename
            + " (" + to_string(num_props) + ")");
      }
      Term prop = propvec[prop_idx];
      Property p(rts, prop);
      std::vector<UnorderedTermMap> cex;
      r = check_prop(engine, bound, p, s, second_solver, cex);

      logger.log(0, "Property {} is {}", prop_idx, to_string(r));

      if (r == FALSE) {
        for (size_t t = 0; t < cex.size(); t++) {
          cout << "AT TIME " << t << endl;
          for (auto elem : cex[t]) {
            cout << "\t" << elem.first << " : " << elem.second << endl;
          }
        }
      }
    } else {
      throw CosaException("Unrecognized file extension " + file_ext
                          + " for file " + filename);
    }
  }
  catch (CosaException & ce) {
    cout << ce.what() << endl;
    cout << "unknown" << endl;
    cout << "b" << prop_idx << endl;
  }
  catch (SmtException & se) {
    cout << se.what() << endl;
    cout << "unknown" << endl;
    cout << "b" << prop_idx << endl;
  }
  catch (std::exception & e) {
    cout << "Caught generic exception..." << endl;
    cout << e.what() << endl;
    cout << "unknown" << endl;
    cout << "b" << prop_idx << endl;
  }

  return status_code;
}
