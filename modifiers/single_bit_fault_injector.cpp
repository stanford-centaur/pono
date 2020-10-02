/*********************                                                        */
/*! \file single_bit_fault_injector.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Modifies a transition system such that valid traces include
**        ones where a single bit can be flipped
**
**/

#include "single_bit_fault_injector.h"

using namespace smt;
using namespace std;

namespace pono {

void SingleBitFaultInjector::create_fault_vals()
{
  const SmtSolver & solver = faulty_fts_.solver();
  unordered_map<Sort, Term> bitflip_cache;

  Term bitflipper;
  Term faultval;
  Term st;
  Term at_most_one;
  Sort st_sort;
  for (auto elem : fts_.state_updates()) {
    st = elem.first;
    st_sort = st->get_sort();

    if (statevars_to_ignore_.find(st) != statevars_to_ignore_.end()) {
      continue;
    }

    // TODO: should we allow flipping boolean state vars?
    //       if there's a witness for the property it could just flip that
    if (st_sort->get_sort_kind() != BV) {
      continue;
    }

    if (bitflip_cache.find(st_sort) != bitflip_cache.end()) {
      bitflipper = bitflip_cache.at(st_sort);
    } else {
      // create a new bitflip val
      bitflipper =
          faulty_fts_.make_statevar(st_sort->to_string() + "_bitflip", st_sort);
      // make it frozen -- the fault model only allows one bitflip
      faulty_fts_.assign_next(bitflipper, bitflipper);
      // constrain it to have at most one bit high
      // bithack x & x - 1 == 0 is equivalent to at most one bit high
      at_most_one = solver->make_term(
          BVAnd,
          bitflipper,
          solver->make_term(BVSub, bitflipper, solver->make_term(1, st_sort)));
      at_most_one =
          solver->make_term(Equal, at_most_one, solver->make_term(0, st_sort));
      faulty_fts_.add_invar(at_most_one);
    }

    // faultval is the original value but with up to one bit flipped
    faultval = solver->make_term(BVXor, st, bitflipper);

    // populate the member variables used in do_injection
    state2faultval_[st] = faultval;
    faultval2state_[faultval] = st;
  }
}

void SingleBitFaultInjector::constrain_to_single_fault()
{
  // at most one for fault sel
  for (auto e1 : faultsel2state_) {
    for (auto e2 : faultsel2state_) {
      Term sel1 = e1.first;
      Term sel2 = e2.first;
      if (sel1 == sel2) {
        // skip when it's the same selector
        continue;
      }

      // Clause: -sel1 \/ -sel2
      faulty_fts_.add_invar(
          faulty_fts_.make_term(Or,
                                faulty_fts_.make_term(Not, sel1),
                                faulty_fts_.make_term(Not, sel2)));
    }
  }

  // once a fault has occured, don't allow it again
  // NOTE: this could be optimized to only use one faultsig
  //       however, we're sticking with the more general version which could
  //       allow multiple faults when/if there are performance issues we can
  //       have a specialized version we could just change the way faultsigs are
  //       created would require refactoring so this is separate from
  //       do_fault_injection which is implemented in FaultInjection

  // create an Or of all fault sigs to determine if any fault has occurred in
  // the past
  Term fault_in_past = faulty_fts_.make_term(false);
  for (auto e : faultsig2state_) {
    Term sig = e.first;
    fault_in_past = faulty_fts_.make_term(Or, fault_in_past, sig);
  }

  // if a fault has happened in the past, disable fault selectors now
  for (auto e : faultsel2state_) {
    Term sel = e.first;
    faulty_fts_.constrain_inputs(faulty_fts_.make_term(
        Implies, fault_in_past, faulty_fts_.make_term(Not, sel)));
  }
}

}  // namespace pono
