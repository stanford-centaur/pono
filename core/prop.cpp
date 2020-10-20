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

Property::Property(TransitionSystem & ts, const Term & p, std::string name)
  : ts_(ts), prop_(p), name_(name)
{
  // find a name if it wasn't provided
  // if no name is associated with it in the ts, then it will just
  // be the to_string of the term
  if (name_.empty())
  {
    name_ = ts_.get_name(prop_);
  }
  initialize();
}

Property::Property(Property & prop, TermTranslator & tt)
    : ts_(prop.ts_, tt),
      // only need to transfer if solvers are different
      // ts_ constructor does the same thing internally
      prop_((prop.transition_system().solver() == tt.get_solver())
                ? prop.prop_
                : tt.transfer_term(prop.prop_)),
      name_(prop.name_)
{
}

Property::~Property() {}

void Property::initialize()
{
  if (ts_.only_curr(prop_)) {
    // nothing to do
    return;
  }

  // otherwise, need to make sure prop_ is over state vars
  logger.log(1,
             "Got next state or input variables in property. Generating a "
             "monitor state.");

  // TODO: process the property
  Term monitor;
  size_t id = 0;
  while (true) {
    try {
      monitor = ts_.make_statevar("_monitor_" + std::to_string(id),
                                  ts_.make_sort(BOOL));
      break;
    }
    catch (SmtException & e) {
      ++id;
    }
  }

  // monitor starts true
  ts_.constrain_init(monitor);

  if (ts_.no_next(prop_)) {
    ts_.assign_next(monitor, prop_);
    ;
  } else if (!ts_.is_functional()) {
    RelationalTransitionSystem & rts =
        static_cast<RelationalTransitionSystem &>(ts_);
    rts.constrain_trans(rts.make_term(Equal, rts.next(monitor), prop_));
  } else {
    assert(ts_.is_functional());
    throw PonoException(
        "Cannot use next in property of a functional transition system.");
  }
}

}  // namespace pono
