/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Makai Mann
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

#pragma once

#include <string>
#include <unordered_map>

namespace pono {

// Engine options
enum Engine
{
  BMC = 0,
  BMC_SP,
  KIND,
  INTERP
};

static const std::unordered_map<std::string, Engine> str2engine(
    { { "bmc", BMC },
      { "bmc-sp", BMC_SP },
      { "ind", KIND },
      { "interp", INTERP } });

const Engine to_engine(std::string s)
{
  if (str2engine.find(s) != str2engine.end()) {
    return str2engine.at(s);
  } else {
    throw PonoException("Unrecognized engine: " + s);
  }
}

/************************************ Default Values
 * *********************************/
static const Engine default_engine = BMC;
static const unsigned int default_prop_idx = 0;
static const unsigned int default_bound = 10;
static const unsigned int default_verbosity = 0;
/********************************* End Default Values
 * ********************************/
}  // namespace pono
