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

#include "core/prop.h"

#include "assert.h"
#include "core/rts.h"
#include "smt-switch/utils.h"
#include "utils/exceptions.h"
#include "utils/logger.h"

using namespace smt;

namespace pono {

Property::Property(const SmtSolver & s, const Term & p,
                   std::string name)
  : solver_(s), prop_(p), name_(name)
{
}

Property::~Property() {}

}  // namespace pono
