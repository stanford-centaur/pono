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

#pragma once

#include <string>

namespace pono {

typedef enum
{
  PROPERTY_TRUE = 0,
  PROPERTY_FALSE = 1,
  PROPERTY_UNKNOWN = 2,
  ERROR = 3
} PonoResult;

std::string to_string(PonoResult r);

}  // namespace pono
