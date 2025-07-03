/*********************                                                        */
/*! \file control_signals.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Useful functions for adding control signal semantics to a system.
**        For example, toggling a clock or setting a reset sequence.
**
**/

#include "modifiers/control_signals.h"

#include <cassert>
#include <cmath>

using namespace smt;
using namespace std;

namespace pono {

void toggle_clock(TransitionSystem & ts, const Term & clock_symbol)
{
  const SmtSolver & s = ts.solver();
  Sort sort = clock_symbol->get_sort();
  SortKind sk = sort->get_sort_kind();
  Sort one_bit_sort = s->make_sort(BV, 1);

  if (clock_symbol->get_sort() != one_bit_sort && sk != BOOL) {
    throw PonoException("Expecting a boolean or one-bit clock sort.");
  }

  Term zero = s->make_term(0, one_bit_sort);
  Term clk_state = clock_symbol;
  assert(!ts.is_next_var(clk_state));
  if (!ts.is_curr_var(clk_state)) {
    Term clk_state = ts.make_statevar(clock_symbol->to_string() + "__state__",
                                      clock_symbol->get_sort());
    ts.constrain_inputs(s->make_term(Equal, clock_symbol, clk_state));
  }
  assert(sk == BV || sk == BOOL);
  if (sk == BV) {
    ts.constrain_init(s->make_term(Equal, clk_state, zero));
    ts.assign_next(clk_state, s->make_term(BVNot, clk_state));
  } else if (sk == BOOL) {
    ts.constrain_init(s->make_term(Not, clk_state));
    ts.assign_next(clk_state, s->make_term(Not, clk_state));
  }
}

Term add_reset_seq(TransitionSystem & ts,
                   const Term & reset_symbol,
                   unsigned long reset_bnd)
{
  const SmtSolver & s = ts.solver();
  Sort reset_sort = reset_symbol->get_sort();
  SortKind sk = reset_sort->get_sort_kind();

  Sort one_bit_sort = s->make_sort(BV, 1);
  if (sk != BOOL && reset_sort != one_bit_sort) {
    throw PonoException("Unexpected reset symbol sort: "
                        + reset_symbol->get_sort()->to_string());
  }

  uint32_t num_bits = ceil(log2(reset_bnd)) + 1;
  Sort bvsort = s->make_sort(BV, num_bits);
  Term reset_bnd_term = s->make_term(reset_bnd, bvsort);
  Term reset_counter = ts.make_statevar("__internal_cosa2_reset_cnt__", bvsort);

  Term in_reset = s->make_term(BVUlt, reset_counter, reset_bnd_term);
  Term reset_done = s->make_term(Not, in_reset);

  ts.constrain_init(
      s->make_term(Equal, reset_counter, s->make_term(0, bvsort)));
  ts.assign_next(
      reset_counter,
      s->make_term(
          Ite,
          reset_done,
          reset_counter,
          s->make_term(BVAdd, reset_counter, s->make_term(1, bvsort))));

  assert(sk == BOOL || sk == BV);
  Term active_reset = reset_symbol;  // if it's boolean
  if (sk == BV) {
    active_reset =
        s->make_term(Equal, reset_symbol, s->make_term(1, one_bit_sort));
  }

  Term inactive_reset = s->make_term(Not, active_reset);
  ts.constrain_inputs(s->make_term(Implies, in_reset, active_reset));
  ts.constrain_inputs(s->make_term(Implies, reset_done, inactive_reset));

  return reset_done;
}

}  // namespace pono
