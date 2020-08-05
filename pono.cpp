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

#include "smt-switch/boolector_factory.h"
#ifdef WITH_MSAT
#include "smt-switch/msat_factory.h"
#endif

#include "bmc.h"
#include "bmc_simplepath.h"
#include "core/fts.h"
#include "frontends/btor2_encoder.h"
#include "frontends/smv_encoder.h"
#include "interpolantmc.h"
#include "kinduction.h"
#include "mbic3.h"
#include "options/options.h"
#include "printers/btor2_witness_printer.h"
#include "printers/vcd_witness_printer.h"
#include "prop.h"
#include "utils/logger.h"

using namespace pono;
using namespace smt;
using namespace std;


ProverResult check_prop(PonoOptions pono_options,
                        Property & p,
                        SmtSolver & s,
                        SmtSolver & second_solver,
                        std::vector<UnorderedTermMap> & cex)
{
  logger.log(1, "Solving property: {}", p.prop());

  logger.log(3, "INIT:\n{}", p.transition_system().init());
  logger.log(3, "TRANS:\n{}", p.transition_system().trans());

  std::shared_ptr<Prover> prover;
  if (pono_options.engine_ == BMC) {
    prover = std::make_shared<Bmc>(pono_options, p, s);
  } else if (pono_options.engine_ == BMC_SP) {
    prover = std::make_shared<BmcSimplePath>(pono_options, p, s);
  } else if (pono_options.engine_ == KIND) {
    prover = std::make_shared<KInduction>(pono_options, p, s);
  } else if (pono_options.engine_ == INTERP) {
    assert(second_solver != NULL);
    prover = std::make_shared<InterpolantMC>(pono_options, p, s, second_solver);
  } else if (pono_options.engine_ == MBIC3) {
    prover = std::make_shared<ModelBasedIC3>(pono_options, p, s);
  } else {
    throw PonoException("Unimplemented engine.");
  }

  ProverResult r = prover->check_until(pono_options.bound_);
  if (r == FALSE && !pono_options.no_witness_) {
    prover->witness(cex);
  }
  return r;
}

int main(int argc, char ** argv)
{
  PonoOptions pono_options;
  ProverResult res = pono_options.parse_and_set_options(argc, argv);
  if (res == ERROR) return res;
  // expected result returned by option parsing and setting is
  // 'pono::UNKNOWN', indicating that options were correctly set and
  // 'ERROR' otherwise, e.g. wrong command line options or
  // incompatible options were passed
  assert(res == pono::UNKNOWN);

  // set logger verbosity -- can only be set once
  logger.set_verbosity(pono_options.verbosity_);

  try {
    SmtSolver s;
    SmtSolver second_solver;
    if (pono_options.engine_ == INTERP) {
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
    } else if (pono_options.engine_ == MBIC3) {
#ifdef WITH_MSAT
      s = MsatSolverFactory::create(false);
      s->set_opt("incremental", "true");
      s->set_opt("produce-models", "true");
#else
      throw PonoException("Model-based IC3 depends on MathSAT currently.");
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
    string file_ext = pono_options.filename_.substr(
        pono_options.filename_.find_last_of(".") + 1);
    if (file_ext == "btor2" || file_ext == "btor") {
      logger.log(2, "Parsing BTOR2 file: {}", pono_options.filename_);
      FunctionalTransitionSystem fts(s);
      BTOR2Encoder btor_enc(pono_options.filename_, fts);
      const TermVec & propvec = btor_enc.propvec();
      unsigned int num_props = propvec.size();
      if (pono_options.prop_idx_ >= num_props) {
        throw PonoException(
            "Property index " + to_string(pono_options.prop_idx_)
            + " is greater than the number of properties in file "
            + pono_options.filename_ + " (" + to_string(num_props) + ")");
      }
      Term prop = propvec[pono_options.prop_idx_];
      Property p(fts, prop);
      vector<UnorderedTermMap> cex;
      res = check_prop(pono_options, p, s, second_solver, cex);
      // we assume that a prover never returns 'ERROR'
      assert(res != ERROR);

      // print btor output
      if (res == FALSE) {
        cout << "sat" << endl;
        cout << "b" << pono_options.prop_idx_ << endl;
        assert(!pono_options.no_witness_ || !cex.size());
        if (cex.size()) {
          print_witness_btor(btor_enc, cex);
          if (!pono_options.vcd_name_.empty()) {
            VCDWitnessPrinter vcdprinter(fts, cex);
            vcdprinter.DumpTraceToFile(pono_options.vcd_name_);
          }
        }
      } else if (res == TRUE) {
        cout << "unsat" << endl;
        cout << "b" << pono_options.prop_idx_ << endl;
      } else {
        assert(res == pono::UNKNOWN);
        cout << "unknown" << endl;
        cout << "b" << pono_options.prop_idx_ << endl;
      }

    } else if (file_ext == "smv") {
      logger.log(2, "Parsing SMV file: {}", pono_options.filename_);
      RelationalTransitionSystem rts(s);
      SMVEncoder smv_enc(pono_options.filename_, rts);
      const TermVec & propvec = smv_enc.propvec();
      unsigned int num_props = propvec.size();
      if (pono_options.prop_idx_ >= num_props) {
        throw PonoException(
            "Property index " + to_string(pono_options.prop_idx_)
            + " is greater than the number of properties in file "
            + pono_options.filename_ + " (" + to_string(num_props) + ")");
      }
      Term prop = propvec[pono_options.prop_idx_];
      Property p(rts, prop);
      std::vector<UnorderedTermMap> cex;
      res = check_prop(pono_options, p, s, second_solver, cex);
      // we assume that a prover never returns 'ERROR'
      assert(res != ERROR);

      logger.log(
          0, "Property {} is {}", pono_options.prop_idx_, to_string(res));

      if (res == FALSE) {
        assert(!pono_options.no_witness_ || cex.size() == 0);
        for (size_t t = 0; t < cex.size(); t++) {
          cout << "AT TIME " << t << endl;
          for (auto elem : cex[t]) {
            cout << "\t" << elem.first << " : " << elem.second << endl;
          }
        }
        assert(!pono_options.no_witness_ || pono_options.vcd_name_.empty());
        if (!pono_options.vcd_name_.empty()) {
          VCDWitnessPrinter vcdprinter(rts, cex);
          vcdprinter.DumpTraceToFile(pono_options.vcd_name_);
        }
      } else if (res == TRUE) {
        cout << "unsat" << endl;
      } else {
        assert(res == pono::UNKNOWN);
        cout << "unknown" << endl;
      }
    } else {
      throw PonoException("Unrecognized file extension " + file_ext
                          + " for file " + pono_options.filename_);
    }
  }
  catch (PonoException & ce) {
    cout << ce.what() << endl;
    cout << "unknown" << endl;
    cout << "b" << pono_options.prop_idx_ << endl;
  }
  catch (SmtException & se) {
    cout << se.what() << endl;
    cout << "unknown" << endl;
    cout << "b" << pono_options.prop_idx_ << endl;
  }
  catch (std::exception & e) {
    cout << "Caught generic exception..." << endl;
    cout << e.what() << endl;
    cout << "unknown" << endl;
    cout << "b" << pono_options.prop_idx_ << endl;
  }

  return res;
}
