/*********************                                                        */
/*! \file single_bit_flip_fault.cpp
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

#include "single_bit_flip_fault.h"

using namespace smt;
using namespace std;

namespace pono {

void SingleBitFlipFault::create_fault_vals()
{
  SmtSolver & solver = faulty_fts_.solver();
  unordered_map<Sort, Term> bitflip_cache;

  Term bitflipper;
  Term faultval;
  Term st;
  Term at_most_one;
  Sort st_sort;
  for (auto elem : fts_.state_updates()) {
    st = elem.first;
    st_sort = st->get_sort();

    // TODO: should we allow flipping boolean state vars?
    //       if there's a witness for the property it could just flip that
    if (st_sort->get_sort_kind() != BV) {
      continue;
    }

    if (bitflip_cache.find(st_sort) != bitflip_cache.end()) {
      bitflipper = bitflip_cache.at(st_sort);
    }
    else {
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

}  // namespace pono
