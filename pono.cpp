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

#include <csignal>
#include <iostream>
#include "assert.h"

#ifdef WITH_PROFILING
#include <gperftools/profiler.h>
#endif

#include "smt-switch/boolector_factory.h"
#ifdef WITH_MSAT
#include "smt-switch/msat_factory.h"
#endif

#include "core/fts.h"
// TODO hide these in make_prover
#include "engines/ceg_prophecy_arrays.h"
#include "engines/cegar_values.h"
#include "engines/ic3ia.h"
// end TODO
#include "frontends/btor2_encoder.h"
#include "frontends/smv_encoder.h"
#include "modifiers/control_signals.h"
#include "modifiers/mod_init_prop.h"
#include "modifiers/prop_monitor.h"
#include "modifiers/static_coi.h"
#include "options/options.h"
#include "printers/btor2_witness_printer.h"
#include "printers/vcd_witness_printer.h"
#include "smt/available_solvers.h"
#include "utils/logger.h"
#include "utils/make_provers.h"
#include "utils/ts_analysis.h"

// TEMP do array abstraction directly here
#include "modifiers/array_abstractor.h"

using namespace pono;
using namespace smt;
using namespace std;


ProverResult check_prop(PonoOptions pono_options,
                        Property & p,
                        const TransitionSystem & ts,
                        const SmtSolver & s,
                        std::vector<UnorderedTermMap> & cex)
{
  logger.log(1, "Solving property: {}", p.name());

  logger.log(3, "INIT:\n{}", ts.init());
  logger.log(3, "TRANS:\n{}", ts.trans());

  Engine eng = pono_options.engine_;

  std::shared_ptr<Prover> prover;
  if (pono_options.cegp_abs_vals_) {
    if (eng != IC3IA_ENGINE || !pono_options.ceg_prophecy_arrays_) {
      throw PonoException(
          "--cegp-abs-vals only supported with -e ic3ia --ceg-prophecy-arrays");
    }
    prover = std::make_shared<CegarValues<CegProphecyArrays<IC3IA>>>(
        p, ts, s, pono_options);
  } else if (pono_options.ceg_prophecy_arrays_) {
    prover = make_ceg_proph_prover(eng, p, ts, s, pono_options);
  } else {
    prover = make_prover(eng, p, ts, s, pono_options);
  }
  assert(prover);

  // TODO: handle this in a more elegant way in the future
  //       consider calling prover for CegProphecyArrays (so that underlying
  //       model checker runs prove unbounded) or possibly, have a command line
  //       flag to pick between the two
  ProverResult r;
  if (pono_options.engine_ == MSAT_IC3IA)
  {
    // HACK MSAT_IC3IA does not support check_until
    r = prover->prove();
  }
  else
  {
    r = prover->check_until(pono_options.bound_);
  }

  if (r == FALSE && !pono_options.no_witness_) {
    bool success = prover->witness(cex);
    if (!success) {
      logger.log(
          0,
          "Only got a partial witness from engine. Not suitable for printing.");
    }
  } else if (r == TRUE && pono_options.check_invar_) {
    try {
      Term invar = prover->invar();
      bool invar_passes = check_invar(ts, p.prop(), invar);
      std::cout << "Invariant Check " << (invar_passes ? "PASSED" : "FAILED")
                << std::endl;
      if (!invar_passes) {
        // shouldn't return true if invariant is incorrect
        r = ProverResult::UNKNOWN;
      }
    }
    catch (PonoException & e) {
      std::cout << "Engine " << pono_options.engine_
                << " does not support getting the invariant." << std::endl;
    }
  }
  return r;
}

// Note: signal handlers are registered only when profiling is enabled.
void profiling_sig_handler(int sig)
{
  std::string signame;
  switch (sig) {
    case SIGINT: signame = "SIGINT"; break;
    case SIGTERM: signame = "SIGTERM"; break;
    case SIGALRM: signame = "SIGALRM"; break;
    default:
      throw PonoException(
          "Caught unexpected signal"
          "in profiling signal handler.");
  }
  logger.log(0, "\n Signal {} received\n", signame);
#ifdef WITH_PROFILING
  ProfilerFlush();
  ProfilerStop();
#endif
  // Switch back to default handling for signal 'sig' and raise it.
  signal(sig, SIG_DFL);
  raise(sig);
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

  // For profiling: set signal handlers for common signals to abort
  // program.  This is necessary to gracefully stop profiling when,
  // e.g., an external time limit is enforced to stop the program.
  if (!pono_options.profiling_log_filename_.empty()) {
    signal(SIGINT, profiling_sig_handler);
    signal(SIGTERM, profiling_sig_handler);
    signal(SIGALRM, profiling_sig_handler);
#ifdef WITH_PROFILING
    logger.log(
        0, "Profiling log filename: {}", pono_options.profiling_log_filename_);
    ProfilerStart(pono_options.profiling_log_filename_.c_str());
#endif
  }

  try {
    // no logging by default
    // could create an option for logging solvers in the future
    SmtSolver s = create_solver_for(
        pono_options.smt_solver_, pono_options.engine_, false);

    // limitations with COI
    if (pono_options.static_coi_) {
      if (!pono_options.no_witness_) {
        logger.log(
            0,
            "Warning: disabling witness production. Temporary restriction -- "
            "Cannot produce witness with option --static-coi");
        pono_options.no_witness_ = true;
      }
      if (pono_options.mod_init_prop_) {
        // Issue explained here:
        // https://github.com/upscale-project/pono/pull/160 will be resolved
        // once state variables removed by COI are removed from init then should
        // do static-coi BEFORE mod-init-prop
        logger.log(0,
                   "Warning: --mod-init-prop and --static-coi don't work "
                   "well together currently.");
      }
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
      // get property name before it is rewritten
      const string prop_name = fts.get_name(prop);

      if (!pono_options.clock_name_.empty()) {
        Term clock_symbol = fts.lookup(pono_options.clock_name_);
        toggle_clock(fts, clock_symbol);
      }
      if (!pono_options.reset_name_.empty()) {
        std::string reset_name = pono_options.reset_name_;
        bool negative_reset = false;
        if (reset_name.at(0) == '~') {
          reset_name = reset_name.substr(1, reset_name.length() - 1);
          negative_reset = true;
        }
        Term reset_symbol = fts.lookup(reset_name);
        if (negative_reset) {
          SortKind sk = reset_symbol->get_sort()->get_sort_kind();
          reset_symbol = (sk == BV) ? s->make_term(BVNot, reset_symbol)
                                    : s->make_term(Not, reset_symbol);
        }
        Term reset_done =
            add_reset_seq(fts, reset_symbol, pono_options.reset_bnd_);
        // guard the property with reset_done
        prop = fts.solver()->make_term(Implies, reset_done, prop);
      }

      if (pono_options.mod_init_prop_) {
        prop = modify_init_and_prop(fts, prop);
      }

      if (pono_options.static_coi_) {
        /* Compute the set of state/input variables related to the
           bad-state property. Based on that information, rebuild the
           transition relation of the transition system. */
        StaticConeOfInfluence coi(fts, { prop }, pono_options.verbosity_);
      }

      if (!fts.only_curr(prop)) {
        logger.log(1, "Got next state or input variables in property. "
                   "Generating a monitor state.");
        prop = add_prop_monitor(fts, prop);
      }

      if (pono_options.assume_prop_) {
        // NOTE: crucial that mod_init_prop and add_prop_monitor passes are
        // before this pass can't assume the non-delayed prop and also
        // delay it
        prop_in_trans(fts, prop);
      }

      vector<UnorderedTermMap> cex;
      Property p(s, prop, prop_name);
      res = check_prop(pono_options, p, fts, s, cex);
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
            vcdprinter.dump_trace_to_file(pono_options.vcd_name_);
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
      // get property name before it is rewritten
      const string prop_name = rts.get_name(prop);

      if (!pono_options.clock_name_.empty()) {
        Term clock_symbol = rts.lookup(pono_options.clock_name_);
        toggle_clock(rts, clock_symbol);
      }
      if (!pono_options.reset_name_.empty()) {
        Term reset_symbol = rts.lookup(pono_options.reset_name_);
        Term reset_done =
            add_reset_seq(rts, reset_symbol, pono_options.reset_bnd_);
        // guard the property with reset_done
        prop = rts.solver()->make_term(Implies, reset_done, prop);
      }

      if (pono_options.mod_init_prop_) {
        prop = modify_init_and_prop(rts, prop);
      }

      if (pono_options.static_coi_) {
        // NOTE: currently only supports FunctionalTransitionSystem
        // but let StaticConeOfInfluence throw the exception
        // and this will change in the future
        /* Compute the set of state/input variables related to the
           bad-state property. Based on that information, rebuild the
           transition relation of the transition system. */
        StaticConeOfInfluence coi(rts, { prop }, pono_options.verbosity_);
      }

      if (!rts.only_curr(prop)) {
        logger.log(1, "Got next state or input variables in property. "
                   "Generating a monitor state.");
        prop = add_prop_monitor(rts, prop);
      }

      if (pono_options.assume_prop_) {
        // NOTE: crucial that mod_init_prop and add_prop_monitor passes are
        // before this pass can't assume the non-delayed prop and also
        // delay it
        prop_in_trans(rts, prop);
      }

      Property p(s, prop, prop_name);
      std::vector<UnorderedTermMap> cex;
      res = check_prop(pono_options, p, rts, s, cex);
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
          vcdprinter.dump_trace_to_file(pono_options.vcd_name_);
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
    cout << "error" << endl;
    cout << "b" << pono_options.prop_idx_ << endl;
    res = ProverResult::ERROR;
  }
  catch (SmtException & se) {
    cout << se.what() << endl;
    cout << "error" << endl;
    cout << "b" << pono_options.prop_idx_ << endl;
    res = ProverResult::ERROR;
  }
  catch (std::exception & e) {
    cout << "Caught generic exception..." << endl;
    cout << e.what() << endl;
    cout << "error" << endl;
    cout << "b" << pono_options.prop_idx_ << endl;
    res = ProverResult::ERROR;
  }

  if (!pono_options.profiling_log_filename_.empty()) {
#ifdef WITH_PROFILING
    ProfilerFlush();
    ProfilerStop();
#endif
  }

  return res;
}
