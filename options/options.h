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
#include "core/proverresult.h"

namespace pono {

// Engine options
enum Engine
{
  BMC = 0,
  BMC_SP,
  KIND,
  INTERP,
  MBIC3
};

const std::unordered_map<std::string, Engine> str2engine({ { "bmc", BMC },
                                                           { "bmc-sp", BMC_SP },
                                                           { "ind", KIND },
                                                           { "interp", INTERP },
                                                           { "mbic3",
                                                             MBIC3 } });

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
        random_seed_(default_random_seed),
        no_witness_(default_no_witness_),
        smt_solver_(default_smt_solver_),
        ic3_cexgen_(default_ic3_cexgen_),
        ic3_indgen_(default_ic3_indgen_),
        ic3_gen_max_iter_(default_ic3_gen_max_iter_),
        ic3_indgen_mode_(default_ic3_indgen_mode_),
        ic3_functional_preimage_(default_ic3_functional_preimage_)
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
  std::string filename_;
  std::string smt_solver_; ///< underlying smt solver

  // ic3 options
  bool ic3_cexgen_;  ///< generalize counterexamples in IC3
  bool ic3_indgen_;  ///< inductive generalization in IC3
  unsigned int ic3_gen_max_iter_; ///< max iterations in ic3 generalization. 0
                                  ///means unbounded
  unsigned int ic3_indgen_mode_; ///< inductive generalization mode [0,1]
  bool ic3_functional_preimage_; ///< functional preimage in IC3

private:
  // Default options
  static const Engine default_engine_ = BMC;
  static const unsigned int default_prop_idx_ = 0;
  static const unsigned int default_bound_ = 10;
  static const unsigned int default_verbosity_ = 0;
  static const unsigned int default_random_seed = 0;
  static const bool default_no_witness_ = false;
  static const std::string default_smt_solver_;
  static const bool default_ic3_cexgen_ = true;
  static const bool default_ic3_indgen_ = true;
  static const unsigned int default_ic3_gen_max_iter_ = 2;
  static const unsigned int default_ic3_indgen_mode_ = 0;
  static const bool default_ic3_functional_preimage_ = false;
};

}  // namespace pono
