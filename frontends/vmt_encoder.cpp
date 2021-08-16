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
**        See https://vmt-lib.fbk.eu/ for more information.
**
**
**/

#include "frontends/vmt_encoder.h"

using namespace smt;
using namespace std;

namespace pono {

VMTEncoder::VMTEncoder(std::string filename, RelationalTransitionSystem & rts)
    : super(rts.get_solver()), filename_(filename), rts_(rts)
{
  set_logic_all();
  int res = parse(filename_);
  assert(!res);  // 0 means success
}

void VMTEncoder::new_symbol(const std::string & name, const smt::Sort & sort)
{
  super::new_symbol(name, sort);
  if (sort->get_sort_kind() != FUNCTION) {
    // treat as an input variable until given :next
    rts_.add_inputvar(lookup_symbol(name));
  }
}

void VMTEncoder::term_attribute(const Term & term,
                                const string & keyword,
                                const string & value)
{
  if (keyword == "next") {
    Term next_var = lookup_symbol(value);
    if (!next_var) {
      // undeclared next var
      // make a new one
      Sort sort = term->get_sort();
      new_symbol(value, sort);
      next_var = lookup_symbol(value);
    }
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
