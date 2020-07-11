/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Ahmed Irfan, Makai Mann, Florian Lonsing
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

#include "ponoresult.h"
#include <cassert>

namespace pono {

std::string to_string(PonoResult r)
{
  if (r == PROPERTY_TRUE) {
    return "PROPERTY TRUE";
  } else if (r == PROPERTY_FALSE) {
    return "PROPERTY FALSE";
  } else if (r == PROPERTY_UNKNOWN) {
    return "PROPERTY UNKNOWN";
  } else {
    assert(r == ERROR);
    return "ERROR";
  }
}

}  // namespace pono
