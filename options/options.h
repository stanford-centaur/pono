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
#include "utils/ponoresult.h"

namespace pono {

// Engine options
enum Engine
{
  BMC = 0,
  BMC_SP,
  KIND,
  INTERP
};

const std::unordered_map<std::string, Engine> str2engine({ { "bmc", BMC },
                                                           { "bmc-sp", BMC_SP },
                                                           { "ind", KIND },
                                                           { "interp",
                                                             INTERP } });

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
        no_witness_(default_no_witness_)
  {
  }

  ~PonoOptions(){};

  PonoResult parse_and_set_options(int argc, char ** argv);

  Engine to_engine(std::string s);

  // Pono options
  Engine engine_;
  unsigned int prop_idx_;
  unsigned int bound_;
  unsigned int verbosity_;
  bool no_witness_;
  std::string vcd_name_;
  std::string filename_;

 private:
  // Default options
  static const Engine default_engine_ = BMC;
  static const unsigned int default_prop_idx_ = 0;
  static const unsigned int default_bound_ = 10;
  static const unsigned int default_verbosity_ = 0;
  static const bool default_no_witness_ = false;
};

}  // namespace pono
