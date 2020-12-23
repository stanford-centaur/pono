/*********************                                                        */
/*! \file syguspdr_ic3formula_helper.cpp
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

#include "utils/sygus_ic3formula_helper.h"


namespace pono {
namespace syntax_analysis {

void GetVariablesFromIC3Formula(const IC3Formula & f, smt::UnorderedTermSet & out) {
  smt::get_free_symbols(f.term, out);
}

std::string PrintModel(const IC3Formula & f) {
  std::string ret;
  for (const auto & c : f.children) {
    auto str = c->to_string();
    auto eqidx = str.find("= ");
    if (eqidx != std::string::npos) {
      str.replace(eqidx, 2, ""); // remove "= ", o
    }
    if (ret.length() == 0)
      ret = str;
    else
      ret = ret + " , " + str;
  }
  return ret;
}

}  // namespace syntax_analysis
}  // namespace pono
