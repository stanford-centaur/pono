/*********************                                                        */
/*! \file vmt_encoder.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Frontend for the Verification Modulo Theories (VMT) format.
**
**
**/

#include "frontends/vmt_encoder.h"

using namespace smt;
using namespace std;

namespace pono {

VMTEncoder::VMTEncoder(std::string filename, RelationalTransitionSystem & rts)
    : super(solver_), filename_(filename), rts_(rts)
{
  solver_ = rts_.solver();
}

void VMTEncoder::term_attribute(const Term & term,
                                const string & keyword,
                                const string & value)
{
  if (keyword == "next") {
    Term next_var = lookup_symbol(value);
    assert(next_var);
    rts_.add_statevar(term, next_var);
  } else if (keyword == "init") {
    rts_.constrain_init(term);
  } else if (keyword == "trans") {
    rts_.constrain_trans(term);
  } else if (keyword == "invar-property") {
    propvec_.push_back(term);
  } else {
    throw PonoException("Unhandled VMT attribute -- :" + keyword + " " + value);
  }
}

}  // namespace pono
