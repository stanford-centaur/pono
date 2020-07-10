/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Makai Mann, Ahmed Irfan
 ** This file is part of the pono project.
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
#include "interpolantmc.h"
#include "kinduction.h"
#include "printers/btor2_witness_printer.h"
#include "printers/vcd_witness_printer.h"
#include "prop.h"
#include "utils/logger.h"
#include "options/options.h"
#include "utils/ponoresult.h"

using namespace pono;
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
  VERBOSITY,
  VCDNAME,
  NOWITNESS
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
  { NOWITNESS, 0, "", "no-witness", Arg::None, "  --no-witness \tDisable printing of witness." },
  { 0, 0, 0, 0, 0, 0 }
};
/*********************************** end Option Handling setup
 * ***************************************/

ProverResult check_prop(Property & p,
                        SmtSolver & s,
                        SmtSolver & second_solver,
                        std::vector<UnorderedTermMap> & cex)
{
  logger.log(1, "Solving property: {}", p.prop());

  logger.log(3, "INIT:\n{}", p.transition_system().init());
  logger.log(3, "TRANS:\n{}", p.transition_system().trans());

  std::shared_ptr<Prover> prover;
  if (pono_options.engine == BMC) {
    prover = std::make_shared<Bmc>(p, s);
  } else if (pono_options.engine == BMC_SP) {
    prover = std::make_shared<BmcSimplePath>(p, s);
  } else if (pono_options.engine == KIND) {
    prover = std::make_shared<KInduction>(p, s);
  } else if (pono_options.engine == INTERP) {
    assert(second_solver != NULL);
    prover = std::make_shared<InterpolantMC>(p, s, second_solver);
  } else {
    throw PonoException("Unimplemented engine.");
  }

  ProverResult r = prover->check_until(pono_options.bound);
  if (r == FALSE && !pono_options.no_witness) {
    prover->witness(cex);
  }
  return r;
}

int main(int argc, char ** argv)
{
  // Set options via the global PonoOptions object 'pono_options'
  // defined in './options/options.h'.
  PonoResult pono_result = pono_options.parse_and_set_options (argc, argv);
  if (pono_result == ERROR)
    return pono_result;
  assert (pono_result == PROPERTY_UNKNOWN);

  // set logger verbosity -- can only be set once
  logger.set_verbosity(pono_options.verbosity);

  try {
    SmtSolver s;
    SmtSolver second_solver;
    if (pono_options.engine == INTERP) {
#ifdef WITH_MSAT
      // need mathsat for interpolant based model checking
      s = MsatSolverFactory::create(false);
      second_solver = MsatSolverFactory::create_interpolating_solver();
#else
      throw PonoException(
          "Interpolation-based model checking requires MathSAT and "
          "this version of pono is built without MathSAT.\nPlease "
          "setup smt-switch with MathSAT and reconfigure using --with-msat.\n"
          "Note: MathSAT has a custom license and you must assume all "
          "responsibility for meeting the license requirements.");
#endif
    } else {
      // boolector is faster but doesn't support interpolants
      s = BoolectorSolverFactory::create(false);
      s->set_opt("produce-models", "true");
      s->set_opt("incremental", "true");
    }

    // TODO: make this less ugly, just need to keep it in scope if using
    //       it would be better to have a generic encoder
    //       and also only create the transition system once
    ProverResult r;
    string file_ext = pono_options.filename.substr(pono_options.filename.find_last_of(".") + 1);
    if (file_ext == "btor2" || file_ext == "btor") {
      logger.log(2, "Parsing BTOR2 file: {}", pono_options.filename);
      FunctionalTransitionSystem fts(s);
      BTOR2Encoder btor_enc(pono_options.filename, fts);
      const TermVec & propvec = btor_enc.propvec();
      unsigned int num_props = propvec.size();
      if (pono_options.prop_idx >= num_props) {
        throw PonoException(
            "Property index " + to_string(pono_options.prop_idx)
            + " is greater than the number of properties in file " + pono_options.filename
            + " (" + to_string(num_props) + ")");
      }
      Term prop = propvec[pono_options.prop_idx];
      Property p(fts, prop);
      vector<UnorderedTermMap> cex;
      r = check_prop(p, s, second_solver, cex);

      // print btor output
      if (r == FALSE) {
        cout << "sat" << endl;
        cout << "b" << pono_options.prop_idx << endl;
        assert (!pono_options.no_witness || !cex.size());
        if (cex.size()) {
          print_witness_btor(btor_enc, cex);
          if (!pono_options.vcd_name.empty()) {
            VCDWitnessPrinter vcdprinter(fts, cex);
            vcdprinter.DumpTraceToFile(pono_options.vcd_name);
          }
        }
        pono_result = PROPERTY_FALSE;
      } else if (r == TRUE) {
        cout << "unsat" << endl;
        cout << "b" << pono_options.prop_idx << endl;
        pono_result = PROPERTY_TRUE;
      } else {
        cout << "unknown" << endl;
        cout << "b" << pono_options.prop_idx << endl;
        pono_result = PROPERTY_UNKNOWN;
      }

    } else if (file_ext == "smv") {
      logger.log(2, "Parsing SMV file: {}", pono_options.filename);
      RelationalTransitionSystem rts(s);
      SMVEncoder smv_enc(pono_options.filename, rts);
      const TermVec & propvec = smv_enc.propvec();
      unsigned int num_props = propvec.size();
      if (pono_options.prop_idx >= num_props) {
        throw PonoException(
            "Property index " + to_string(pono_options.prop_idx)
            + " is greater than the number of properties in file " + pono_options.filename
            + " (" + to_string(num_props) + ")");
      }
      Term prop = propvec[pono_options.prop_idx];
      Property p(rts, prop);
      std::vector<UnorderedTermMap> cex;
      r = check_prop(p, s, second_solver, cex);
      logger.log(0, "Property {} is {}", pono_options.prop_idx, to_string(r));

      if (r == FALSE) {
        pono_result = PROPERTY_FALSE;
        assert (!pono_options.no_witness || cex.size() == 0);
        for (size_t t = 0; t < cex.size(); t++) {
          cout << "AT TIME " << t << endl;
          for (auto elem : cex[t]) {
            cout << "\t" << elem.first << " : " << elem.second << endl;
          }
        }
        assert (!pono_options.no_witness || pono_options.vcd_name.empty());
        if (!pono_options.vcd_name.empty()) {
          VCDWitnessPrinter vcdprinter(rts, cex);
          vcdprinter.DumpTraceToFile(pono_options.vcd_name);
        }
      }
      else if (r == TRUE) {
        cout << "unsat" << endl;
        pono_result = PROPERTY_TRUE;
      }
      else
        {
          cout << "unknown" << endl;
          pono_result = PROPERTY_UNKNOWN;
        }
    } else {
      throw PonoException("Unrecognized file extension " + file_ext
                          + " for file " + pono_options.filename);
    }
  }
  catch (PonoException & ce) {
    cout << ce.what() << endl;
    cout << "unknown" << endl;
    cout << "b" << pono_options.prop_idx << endl;
  }
  catch (SmtException & se) {
    cout << se.what() << endl;
    cout << "unknown" << endl;
    cout << "b" << pono_options.prop_idx << endl;
  }
  catch (std::exception & e) {
    cout << "Caught generic exception..." << endl;
    cout << e.what() << endl;
    cout << "unknown" << endl;
    cout << "b" << pono_options.prop_idx << endl;
  }

  return pono_result;
}
