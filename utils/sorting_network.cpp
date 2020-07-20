/*********************                                                        */
/*! \file sorting_network.cpp
** \verbatim
** Top contributors (to current version):
**   Ahmed Irfan, Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Sorting networks for symbolically sorting terms.
**
**
**/

#include <cassert>
#include <cstdlib>

#include "sorting_network.h"

using namespace std;
using namespace smt;

namespace cosa {

void BoolSortingNetwork::do_sorting()
{
  if (!boolvec_.size()) {
    return;
  }
  sorted_boolvec_ = sorting_network_helper(boolvec_);
}

TermVec BoolSortingNetwork::sorting_network_helper(TermVec boolvec) const
{
  if (boolvec.size() <= 1) {
    return boolvec;
  } else if (boolvec.size() == 2) {
    return two_comparator(boolvec[0], boolvec[1]);
  } else {
    size_t pivot = boolvec.size() / 2;
    // TODO: Check this
    TermVec left_boolvec = sorting_network_helper(
        TermVec(boolvec.begin(), boolvec.begin() + pivot));
    TermVec right_boolvec =
        sorting_network_helper(TermVec(boolvec.begin() + pivot, boolvec.end()));
    return merge(left_boolvec, right_boolvec);
  }
}

TermVec BoolSortingNetwork::two_comparator(Term bool0, Term bool1) const
{
  return TermVec({ solver_->make_term(Or, bool0, bool1),
                   solver_->make_term(And, bool0, bool1) });
}

TermVec BoolSortingNetwork::merge(TermVec boolvec0, TermVec boolvec1) const
{
  size_t boolvec0_size = boolvec0.size();
  size_t boolvec1_size = boolvec1.size();

  // TODO check this -- probably not true
  // assert(boolvec0_size < boolvec1_size);
  // TODO this one probably is true
  assert(abs(((int)boolvec0_size) - ((int)boolvec1_size)) <= 1);

  // base cases
  if (!boolvec0_size) {
    return boolvec1;
  } else if (!boolvec1_size) {
    return boolvec0;
  } else if (boolvec0_size == 1 && boolvec1_size == 1) {
    return two_comparator(boolvec0[0], boolvec1[0]);
  }

  // main, recursive merging logic
  bool boolvec0_is_even = (boolvec0_size % 2) == 0;
  bool boolvec1_is_even = (boolvec1_size % 2) == 0;

  // normalize order
  if (!boolvec0_is_even && boolvec1_is_even) {
    return merge(boolvec1, boolvec0);
  }

  // assert(boolvec0_size <= boolvec1_size);

  TermVec boolvec0_odd;
  TermVec boolvec0_even;

  for (size_t i = 0; i < boolvec0.size(); i++) {
    Term elem = boolvec0[i];
    if ((i + 1) % 2 == 0) {
      boolvec0_even.push_back(elem);
    } else {
      boolvec0_odd.push_back(elem);
    }
  }

  TermVec boolvec1_odd;
  TermVec boolvec1_even;

  for (size_t i = 0; i < boolvec1.size(); i++) {
    Term elem = boolvec1[i];
    if ((i + 1) % 2 == 0) {
      boolvec1_even.push_back(elem);
    } else {
      boolvec1_odd.push_back(elem);
    }
  }

  TermVec merged_odd = merge(boolvec0_odd, boolvec1_odd);
  TermVec merged_even = merge(boolvec0_even, boolvec1_even);

  size_t merged_odd_size = merged_odd.size();
  size_t merged_even_size = merged_even.size();

  TermVec output;
  // cases
  // even and even
  if (boolvec0_is_even && boolvec1_is_even) {
    output.push_back(merged_odd[0]);

    if (merged_odd_size > 0) {
      for (size_t i = 0; i < merged_odd_size - 1; i++) {
        for (Term el : two_comparator(merged_odd[i + 1], merged_even[i])) {
          output.push_back(el);
        }
      }
    }
    output.push_back(merged_even[merged_even_size - 1]);
  }
  // even and odd
  else if (boolvec0_is_even && !boolvec1_is_even) {
    output.push_back(merged_odd[0]);
    for (size_t i = 0; i < merged_even_size; i++) {
      for (Term el : two_comparator(merged_odd[i + 1], merged_even[i])) {
        output.push_back(el);
      }
    }
  }
  // odd and odd
  else if (!boolvec0_is_even && !boolvec1_is_even) {
    output.push_back(merged_odd[0]);
    for (size_t i = 0; i < merged_even_size; i++) {
      for (Term el : two_comparator(merged_odd[i + 1], merged_even[i])) {
        output.push_back(el);
      }
    }
    output.push_back(merged_odd[merged_odd_size - 1]);
  } else {
    // this case should not occur (normalized order)
    assert(false);
  }

  assert(boolvec0_size + boolvec1_size == output.size());

  return output;
}

}  // namespace cosa
