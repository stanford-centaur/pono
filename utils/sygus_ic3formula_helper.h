/*********************                                                        */
/*! \file syguspdr_ic3formula_helper.h
** \verbatim
** Top contributors (to current version):
**   Hongce Zhang
** This file is part of the pono project.
** Copyright (c) 2020 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Utility functions to get variables from ic3formula
**        and to print it
**
**/

#pragma once

#include "engines/ic3base.h"

#include <string>

namespace pono {

namespace syntax_analysis {

void GetVariablesFromIC3Formula(const IC3Formula & f, smt::UnorderedTermSet & out);

std::string PrintModel(const IC3Formula & f);

}  // namespace syntax_analysis
}  // namespace pono


