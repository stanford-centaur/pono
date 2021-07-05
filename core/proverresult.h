/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Ahmed Irfan, Makai Mann
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

#include <iostream>
#include <string>

namespace pono {

typedef enum
{
  UNKNOWN = -1,
  FALSE = 0,
  TRUE = 1,
  ERROR = 2
} ProverResult;

std::string to_string(ProverResult r);

std::ostream & operator<<(std::ostream & o, ProverResult r);

}  // namespace pono
