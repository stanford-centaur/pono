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

#include <iostream>
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
  IC3_BITS,
  MBIC3,
  IC3IA_ENGINE,
  MSAT_IC3IA,
  IC3SA_ENGINE,
  SYGUS_PDR
  // NOTE: if adding an IC3 variant,
  // make sure to update ic3_variants_set in options/options.cpp
  // used for setting solver options appropriately
};

// NOTE keep this up to date so that setting
//      solver options based on the engine works
//      as expected.
//      IC3 uses the solver differently than other
//      techniques and thus different options are
//      appropriate
/** Returns the set of all IC3 variant engines */
const std::unordered_set<Engine> & ic3_variants();

const std::unordered_map<std::string, Engine> str2engine(
    { { "bmc", BMC },
      { "bmc-sp", BMC_SP },
      { "ind", KIND },
      { "interp", INTERP },
      { "mbic3", MBIC3 },
      { "ic3bool", IC3_BOOL},
      { "ic3bits", IC3_BITS },
      { "ic3ia", IC3IA_ENGINE },
      { "msat-ic3ia", MSAT_IC3IA },
      { "ic3sa", IC3SA_ENGINE },
      { "sygus-pdr", SYGUS_PDR } });

// SyGuS mode option
enum SyGuSTermMode{
  FROM_DESIGN_LEARN_EXT = 0,
  VAR_C_EXT = 1,
  SPLIT_FROM_DESIGN = 2,
  VAR_C_EQ_LT = 3,
  TERM_MODE_AUTO = 4
};

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
        witness_(default_witness_),
        reset_bnd_(default_reset_bnd_),
        random_seed_(default_random_seed),
        smt_solver_(default_smt_solver_),
        logging_smt_solver_(default_logging_smt_solver_),
        static_coi_(default_static_coi_),
        show_invar_(default_show_invar_),
        check_invar_(default_check_invar_),
        ic3_pregen_(default_ic3_pregen_),
        ic3_indgen_(default_ic3_indgen_),
        ic3_gen_max_iter_(default_ic3_gen_max_iter_),
        mbic3_indgen_mode(default_mbic3_indgen_mode),
        ic3_functional_preimage_(default_ic3_functional_preimage_),
        ic3_unsatcore_gen_(default_ic3_unsatcore_gen_),
        ic3ia_reduce_preds_(default_ic3ia_reduce_preds_),
        ic3ia_track_important_vars_(default_ic3ia_track_important_vars_),
        ic3sa_func_refine_(default_ic3sa_func_refine_),
        profiling_log_filename_(default_profiling_log_filename_),
        pseudo_init_prop_(default_pseudo_init_prop_),
        assume_prop_(default_assume_prop_),
        ceg_prophecy_arrays_(default_ceg_prophecy_arrays_),
        cegp_timed_axiom_red_(default_cegp_timed_axiom_red_),
        cegp_consec_axiom_red_(default_cegp_consec_axiom_red_),
        cegp_nonconsec_axiom_red_(default_cegp_nonconsec_axiom_red_),
        cegp_force_restart_(default_cegp_force_restart_),
        cegp_abs_vals_(default_cegp_abs_vals_),
        cegp_abs_vals_cutoff_(default_cegp_abs_vals_cutoff_),
        cegp_strong_abstraction_(default_cegp_strong_abstraction_),
        ceg_bv_arith_(default_ceg_bv_arith_),
        ceg_bv_arith_min_bw_(default_ceg_bv_arith_min_bw_),
        promote_inputvars_(default_promote_inputvars_),
        sygus_term_mode_(default_sygus_term_mode_),
        sygus_term_extract_depth_(default_sygus_term_extract_depth_),
        sygus_initial_term_width_(default_sygus_initial_term_width_),
        sygus_initial_term_inc_(default_sygus_initial_term_inc_),
        sygus_accumulated_term_bound_(default_sygus_accumulated_term_bound_),
        sygus_use_operator_abstraction_(
            default_sygus_use_operator_abstraction_),
        ic3sa_initial_terms_lvl_(default_ic3sa_initial_terms_lvl_),
        ic3sa_interp_(default_ic3sa_interp_),
        print_wall_time_(default_print_wall_time_),
        bmc_bound_start_(default_bmc_bound_start_),
        bmc_bound_step_(default_bmc_bound_step_),
        bmc_neg_init_step_(default_bmc_neg_init_step_),
        bmc_exponential_step_(default_bmc_exponential_step_),
        bmc_single_bad_state_(default_bmc_single_bad_state_),
        bmc_neg_bad_step_(default_bmc_neg_bad_step_),
        bmc_neg_bad_step_all_(default_bmc_neg_bad_step_all_),
        bmc_min_cex_linear_search_(default_bmc_min_cex_linear_search_),
        bmc_min_cex_less_inc_bin_search_(default_bmc_min_cex_less_inc_bin_search_),
        bmc_allow_non_minimal_cex_(default_bmc_allow_non_minimal_cex_),
        kind_no_simple_path_check_(default_kind_no_simple_path_check_),
        kind_eager_simple_path_check_(default_kind_eager_simple_path_check_),
        kind_no_multi_call_simple_path_check_(default_kind_no_multi_call_simple_path_check_),
        kind_no_ind_check_init_states_(default_kind_no_ind_check_init_states_),
        kind_no_ind_check_(default_kind_no_ind_check_),
        kind_no_ind_check_property_(default_kind_no_ind_check_property_),
        kind_one_time_base_check_(default_kind_one_time_base_check_),
        kind_bound_step_(default_kind_bound_step_)
  {
  }

  ~PonoOptions(){};

  /** Parse and set options given argc and argv from main
   *  @param argc
   *  @param argv
   *  @param expect_file if true expects a filename to read
   */
  ProverResult parse_and_set_options(int argc,
                                     char ** argv,
                                     bool expect_file = true);

  /** Parse and set options given vector of options
   *  @param opts vector of command line options
   *  @param expect_file if true expects a filename to read
   */
  ProverResult parse_and_set_options(std::vector<std::string> & opts,
                                     bool expect_file = true);

  Engine to_engine(std::string s);

  // Pono options
  Engine engine_;
  unsigned int prop_idx_;
  unsigned int bound_;
  unsigned int verbosity_;
  unsigned int random_seed_;
  bool witness_;
  std::string vcd_name_;
  std::string reset_name_;
  size_t reset_bnd_;
  std::string clock_name_;
  std::string filename_;
  smt::SolverEnum smt_solver_;  ///< underlying smt solver
  bool logging_smt_solver_;
  bool static_coi_;
  bool show_invar_;   ///< display invariant when running from command line
  bool check_invar_;  ///< check invariants (if available) when run through CLI
  // ic3 options
  bool ic3_pregen_;  ///< generalize counterexamples in IC3
  bool ic3_indgen_;  ///< inductive generalization in IC3
  unsigned int ic3_gen_max_iter_; ///< max iterations in ic3 generalization. 0
                                  ///means unbounded
  unsigned int mbic3_indgen_mode;  ///< inductive generalization mode [0,2]
  bool ic3_functional_preimage_; ///< functional preimage in IC3
  bool ic3_unsatcore_gen_;  ///< generalize a cube during relative inductiveness
                            ///< check with unsatcore
  bool ic3ia_reduce_preds_;  ///< reduce predicates with unsatcore in IC3IA
  bool ic3ia_track_important_vars_;  ///< prioritize predicates with marked
                                     ///< important variables
  bool ic3sa_func_refine_;  ///< try functional unrolling in refinement
  std::string profiling_log_filename_;
  bool pseudo_init_prop_;  ///< replace init and prop with boolean state vars
  bool assume_prop_;       ///< assume property in pre-state
  // ceg-prophecy-arrays options
  bool ceg_prophecy_arrays_;
  bool cegp_timed_axiom_red_;      ///< reduce axioms with an unsat core when
                                   ///< enumerating
  bool cegp_consec_axiom_red_;     ///< reduce consecutive axioms before lifting
  bool cegp_nonconsec_axiom_red_;  ///< reduce nonconsecutive axioms before
                                   ///< prophecizing
  bool cegp_force_restart_;  ///< force underlying engine to restart after
                             ///< refinement
  bool cegp_abs_vals_;  ///< abstract values on top of ceg-prophecy-arrays
  size_t cegp_abs_vals_cutoff_;  ///< cutoff to abstract a value
  bool cegp_strong_abstraction_;  ///< use strong abstraction (no equality UFs)
  bool ceg_bv_arith_;            ///< CEGAR -- Abstract BV arithmetic operators
  size_t ceg_bv_arith_min_bw_;   ///< Only abstract operators having bitwidth
                                 ///< strictly greater than this number
  bool promote_inputvars_;
  // sygus-pdr options
  SyGuSTermMode sygus_term_mode_; ///< SyGuS term production mode
  unsigned sygus_term_extract_depth_; ///< SyGuS Term extraction depth for existing terms
  unsigned sygus_initial_term_width_; ///< SyGuS Control and data width seperator
  unsigned sygus_initial_term_inc_; ///< SyGuS Control and data width seperator increment bound
  unsigned sygus_accumulated_term_bound_; ///< SyGuS Term accumulation bound count
  unsigned sygus_use_operator_abstraction_; ///< SyGuS abstract and avoid use some operators
  size_t ic3sa_initial_terms_lvl_;  ///< configures where to find terms for
                                    ///< initial abstraction
  bool ic3sa_interp_;
  // print wall clock time spent in entire execution
  bool print_wall_time_;

  // BMC interval options (these options are modifiers of the 'BMC' engine;
  //   they do not apply to engine 'BMC-SP')
  // Default bmc_bound_start_ == 0, which starts search for cex at
  // unrolling depth 0 like traditional BMC.
  unsigned bmc_bound_start_;
  // Default: bmc_bound_step_ == 1, which results in traditional BMC
  // where every bound is checked one by one. bmc_bound_step_ is the
  // value by which the current unrolling depth is increased. For
  // bmc_bound_step_ > 1, BMC searches for cex in intervals of size
  // bmc_bound_step_.
  unsigned bmc_bound_step_;
  // BMC: add negated initial state predicate in steps k > 0 (default: false)
  bool bmc_neg_init_step_;
  // BMC: double the bound in each step starting from
  // 'bmc_bound_start_'; if bmc_bound_start_ == 0, this results in
  // exploration of bounds 0,1,2,4,8,...
  bool bmc_exponential_step_;
  // BMC EXPERT OPTION: do not add a disjunctive bad state property
  // representing an interval, but a single bad state literal at bound k;
  bool bmc_single_bad_state_;
  // BMC: add negated bad state predicate depending on reached_k_ (default: false)
  bool bmc_neg_bad_step_;
  // BMC: like 'bmc_neg_bad_step_' but adds negated bad state predicate for all
  // seen bounds (default: false)
  bool bmc_neg_bad_step_all_;
  // Apply linear instead of binary search for minimal counterexample
  // after a counterexample was found within an interval
  bool bmc_min_cex_linear_search_;
  // BMC: apply less incremental binary search
  bool bmc_min_cex_less_inc_bin_search_;
  // BMC: when using disjunctive bad state property, skip the
  // minimization phase after an interval containing a cex was found,
  // i.e., skip binary or linear search for shortest cex in that
  // interval
  bool bmc_allow_non_minimal_cex_;
  // K-induction: omit simple path check (might cause incompleteness)
  bool kind_no_simple_path_check_;
  // K-induction: eager simple path check (default: lazy check)
  bool kind_eager_simple_path_check_;
  // K-induction: no multi-call simple path check
  bool kind_no_multi_call_simple_path_check_;
  // K-induction: skip inductive case check based on initial states
  bool kind_no_ind_check_init_states_;
  // K-induction: skip inductive case check (EXPERT OPTION: will cause
  // incompleteness for most problem instances); this option implies
  // 'kind_no_ind_check_init_states_ == true' and
  // 'kind_no_ind_check_property_ == true'
  bool kind_no_ind_check_;
  // K-induction: skip inductive case check based on property (EXPERT
  // OPTION: will cause incompleteness for most problem instances)
  bool kind_no_ind_check_property_;
  // K-induction: check base case only once after inductive case check was unsatisfiable
  bool kind_one_time_base_check_;
  // K-induction: amount of steps by which transition relation is unrolled
  unsigned kind_bound_step_;

private:
  // Default options
  static const Engine default_engine_ = BMC;
  static const unsigned int default_prop_idx_ = 0;
  static const unsigned int default_bound_ = 10;
  static const unsigned int default_verbosity_ = 0;
  static const unsigned int default_random_seed = 0;
  static const bool default_witness_ = false;
  static const bool default_static_coi_ = false;
  static const bool default_show_invar_ = false;
  static const bool default_check_invar_ = false;
  static const size_t default_reset_bnd_ = 1;
  // TODO distinguish when solver is not set and choose a
  //      good solver for the provided engine automatically
  static const smt::SolverEnum default_smt_solver_ = smt::BTOR;
  static const bool default_logging_smt_solver_ = false;
  static const bool default_ic3_pregen_ = true;
  static const bool default_ic3_indgen_ = true;
  static const unsigned int default_ic3_gen_max_iter_ = 2;
  static const unsigned int default_mbic3_indgen_mode = 0;
  static const bool default_ic3_functional_preimage_ = false;
  static const bool default_ic3_unsatcore_gen_ = true;
  static const bool default_ic3ia_reduce_preds_ = true;
  static const bool default_ic3ia_track_important_vars_ = true;
  static const bool default_ic3sa_func_refine_ = true;
  static const std::string default_profiling_log_filename_;
  static const bool default_pseudo_init_prop_ = false;
  static const bool default_assume_prop_ = false;
  static const bool default_ceg_prophecy_arrays_ = false;
  static const bool default_cegp_timed_axiom_red_ = true;
  static const bool default_cegp_consec_axiom_red_ = true;
  static const bool default_cegp_nonconsec_axiom_red_ = true;
  static const bool default_cegp_force_restart_ = false;
  static const bool default_cegp_abs_vals_ = false;
  static const size_t default_cegp_abs_vals_cutoff_ = 100;
  static const bool default_cegp_strong_abstraction_ = false;
  static const bool default_ceg_bv_arith_ = false;
  static const size_t default_ceg_bv_arith_min_bw_ = 16;
  static const bool default_promote_inputvars_ = false;
  static const SyGuSTermMode default_sygus_term_mode_ = TERM_MODE_AUTO;
  static const unsigned default_sygus_term_extract_depth_ = 0;
  static const unsigned default_sygus_initial_term_width_ = 8;
  static const unsigned default_sygus_initial_term_inc_ = 8;
  static const unsigned default_sygus_accumulated_term_bound_ = 0;
  static const unsigned default_sygus_use_operator_abstraction_ = 0;
  // default is the highest level
  static const size_t default_ic3sa_initial_terms_lvl_ = 4;
  static const bool default_ic3sa_interp_ = false;
  static const bool default_print_wall_time_ = false;
  static const unsigned default_bmc_bound_start_ = 0;
  static const unsigned default_bmc_bound_step_ = 1;
  static const bool default_bmc_neg_init_step_ = false;
  static const bool default_bmc_exponential_step_ = false;
  static const bool default_bmc_single_bad_state_ = false;
  static const bool default_bmc_neg_bad_step_ = false;
  static const bool default_bmc_neg_bad_step_all_ = false;
  static const bool default_bmc_min_cex_linear_search_ = false;
  static const bool default_bmc_min_cex_less_inc_bin_search_ = false;
  static const bool default_bmc_allow_non_minimal_cex_ = false;
  static const bool default_kind_no_simple_path_check_ = false;
  static const bool default_kind_eager_simple_path_check_ = false;
  static const bool default_kind_no_multi_call_simple_path_check_ = false;
  static const bool default_kind_no_ind_check_init_states_ = false;
  static const bool default_kind_no_ind_check_ = false;
  static const bool default_kind_no_ind_check_property_ = false;
  static const bool default_kind_one_time_base_check_ = false;
  static const unsigned default_kind_bound_step_ = 1;
};

// Useful functions for printing etc...

std::string to_string(Engine e);

std::ostream & operator<<(std::ostream & o, Engine e);

}  // namespace pono
