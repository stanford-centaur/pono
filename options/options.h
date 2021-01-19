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

#pragma once

#include <unordered_map>
#include <unordered_set>

#include "core/proverresult.h"
#include "smt-switch/smt.h"

namespace pono {

// Engine options
enum Engine
{
  NONE = -1,
  BMC = 0,
  BMC_SP,
  KIND,
  INTERP,
  IC3_BOOL,
  MBIC3,
  IC3IA_ENGINE,
  MSAT_IC3IA
  // NOTE: if adding an IC3 variant,
  // make sure to update ic3_variants in smt/available_solvers.cpp
  // used for setting solver options appropriately
};

const std::unordered_map<std::string, Engine> str2engine(
    { { "bmc", BMC },
      { "bmc-sp", BMC_SP },
      { "ind", KIND },
      { "interp", INTERP },
      { "mbic3", MBIC3 },
      { "ic3ia", IC3IA_ENGINE },
      { "msat-ic3ia", MSAT_IC3IA } });

/*************************************** Options class
 * ************************************************/

class PonoOptions
{
  // TODO: for now, all class members (i.e., Pono options) are public
  // and set via a call of 'parse_and_set_options()'. For more
  // flexibility, add get/set functions.

 public:
  PonoOptions()
      : engine_(default_engine_),
        prop_idx_(default_prop_idx_),
        bound_(default_bound_),
        verbosity_(default_verbosity_),
        no_witness_(default_no_witness_),
        reset_bnd_(default_reset_bnd_),
        random_seed_(default_random_seed),
        smt_solver_(default_smt_solver_),
        static_coi_(default_static_coi_),
        check_invar_(default_check_invar_),
        ic3_pregen_(default_ic3_pregen_),
        ic3_indgen_(default_ic3_indgen_),
        ic3_gen_max_iter_(default_ic3_gen_max_iter_),
        ic3_reset_interval_(default_ic3_reset_interval_),
        mbic3_indgen_mode(default_mbic3_indgen_mode),
        ic3_functional_preimage_(default_ic3_functional_preimage_),
        ceg_prophecy_arrays_(default_ceg_prophecy_arrays_),
        cegp_axiom_red_(default_cegp_axiom_red_),
        profiling_log_filename_(default_profiling_log_filename_),
        mod_init_prop_(default_mod_init_prop_),
        assume_prop_(default_assume_prop_),
        cegp_abs_vals_(default_cegp_abs_vals_),
        cegp_abs_vals_cutoff_(default_cegp_abs_vals_cutoff_)
  {
  }

  ~PonoOptions(){};

  ProverResult parse_and_set_options(int argc, char ** argv);

  Engine to_engine(std::string s);

  // Pono options
  Engine engine_;
  unsigned int prop_idx_;
  unsigned int bound_;
  unsigned int verbosity_;
  unsigned int random_seed_;
  bool no_witness_;
  std::string vcd_name_;
  std::string reset_name_;
  size_t reset_bnd_;
  std::string clock_name_;
  std::string filename_;
  smt::SolverEnum smt_solver_;  ///< underlying smt solver
  bool static_coi_;
  bool check_invar_;  ///< check invariants (if available) when run through CLI
  // ic3 options
  bool ic3_pregen_;  ///< generalize counterexamples in IC3
  bool ic3_indgen_;  ///< inductive generalization in IC3
  unsigned int ic3_reset_interval_;  ///< number of check sat calls before
                                     ///< resetting. 0 means unbounded
  unsigned int ic3_gen_max_iter_; ///< max iterations in ic3 generalization. 0
                                  ///means unbounded
  unsigned int mbic3_indgen_mode;  ///< inductive generalization mode [0,2]
  bool ic3_functional_preimage_; ///< functional preimage in IC3
  // ceg-prophecy-arrays options
  bool ceg_prophecy_arrays_;
  bool cegp_axiom_red_;  ///< reduce axioms with an unsat core in ceg prophecy
  std::string profiling_log_filename_;
  bool mod_init_prop_;  ///< replace init and prop with boolean state vars
  bool assume_prop_;    ///< assume property in pre-state
  bool cegp_abs_vals_;  ///< abstract values on top of ceg-prophecy-arrays
  size_t cegp_abs_vals_cutoff_;  ///< cutoff to abstract a value

 private:
  // Default options
  static const Engine default_engine_ = BMC;
  static const unsigned int default_prop_idx_ = 0;
  static const unsigned int default_bound_ = 10;
  static const unsigned int default_verbosity_ = 0;
  static const unsigned int default_random_seed = 0;
  static const bool default_no_witness_ = false;
  static const bool default_ceg_prophecy_arrays_ = false;
  static const bool default_static_coi_ = false;
  static const bool default_check_invar_ = false;
  static const size_t default_reset_bnd_ = 1;
  // TODO distinguish when solver is not set and choose a
  //      good solver for the provided engine automatically
  static const smt::SolverEnum default_smt_solver_ = smt::BTOR;
  static const bool default_ic3_pregen_ = true;
  static const bool default_ic3_indgen_ = true;
  static const unsigned int default_ic3_reset_interval_ = 5000;
  static const unsigned int default_ic3_gen_max_iter_ = 2;
  static const unsigned int default_mbic3_indgen_mode = 0;
  static const bool default_ic3_functional_preimage_ = false;
  static const bool default_cegp_axiom_red_ = true;
  static const std::string default_profiling_log_filename_;
  static const bool default_mod_init_prop_ = false;
  static const bool default_assume_prop_ = true;
  static const bool default_cegp_abs_vals_ = false;
  static const size_t default_cegp_abs_vals_cutoff_ = 100;
};

}  // namespace pono
