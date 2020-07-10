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

#include "utils/ponoresult.h"
#include "utils/exceptions.h"
#include "options/defaults.h"

/*************************************** Options class
 * ************************************************/
// Meant to be used as a singleton class

namespace pono {

class PonoOptions
{

  // TODO: for now, all class members (i.e., Pono options) are public
  // and set via a call of 'parse_and_set_options()'. For more
  // flexibility, add get/set functions.
  
 public:

 PonoOptions() :
  engine(default_engine),
    prop_idx(default_prop_idx),
    bound(default_bound),
    verbosity(default_verbosity),
    no_witness(default_no_witness)
    {
      if (!instance_created)
        instance_created = true;
      else
        throw PonoException("Cannot create objects of class PonoOptions, "\
                            "use global object 'pono_options' instead.");
    }

  ~PonoOptions() {};
  
  PonoResult parse_and_set_options (int argc, char ** argv);

  // Pono options
  Engine engine;
  unsigned int prop_idx;
  unsigned int bound;
  unsigned int verbosity;
  bool no_witness;
  std::string vcd_name;
  std::string filename;
  
 private:
  // flag to prevent creation of more than one class object
  static bool instance_created;
};

// globally available options instance
extern PonoOptions pono_options;

}  // namespace pono
