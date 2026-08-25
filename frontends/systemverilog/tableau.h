/*!
 * \file tableau.h
 * \brief Pure tableau/latch-building primitives for SVA/LTL encoding.
 * \author Áron Ricardo Perez-Lopez
 * \date 2026
 * \copyright See the LICENSE file in the top-level source directory.
 *
 * Tableau holds the SVA/LTL encoder's tableau/latch-building gadgets:
 * history-chain latching for `$past`/delayed-boolean tracking
 * (make_history_chain, delay_bool), the "first cycle" and "before cycle k"
 * flags (init_flag, before_cycle), the `disable iff` shift-window OR
 * (disable_window), and the one-step "promise" testers for the LTL tableau's
 * X/G/F/R/U operators (make_X/G/F/R/U), with justice (fairness) conditions
 * for the least-fixpoint operators F and (strong) U.
 *
 * Every method here only needs the transition system and solver, plus its
 * own latch-naming counter and cached "first cycle" flag -- none of them
 * ever walk a slang AST or convert a slang Expression to a Term. That work
 * (assertion_expr_to_bool(), offsets_ending_now(), ltl_to_sat(), and
 * friends -- see sva.cpp) stays on SystemVerilogEncoder, since it needs
 * expr_encoder_.expr_to_term(), which in turn needs SymbolTable's
 * declaration/wire-resolution state. Tableau is instead used by
 * SystemVerilogEncoder as an ordinary owned member: it holds no reference
 * back to SystemVerilogEncoder, and takes no callbacks -- every call passes
 * the Terms/counts/name prefix it needs explicitly.
 */
#pragma once

#include <cstdint>
#include <string>

#include "core/fts.h"
#include "smt-switch/smt.h"

namespace pono {

class Tableau
{
 public:
  Tableau(FunctionalTransitionSystem & fts, const smt::SmtSolver & solver);

  /** Build a chain of `n` 1-cycle latch state vars that track `value`
   *  over n cycles, returning a Term equal to the value from `n` cycles
   *  ago.  Each latch is initialised to 0 and its next-state is the
   *  previous link in the chain.  Used by both `$past(...)` and
   *  sequence-delay assertion expressions.
   *  @param value the current-cycle value to delay
   *  @param n     the number of cycles of delay (n == 0 is a no-op,
   *         returning `value` unchanged)
   *  @param name_prefix the hierarchical name prefix of the module this
   *         latch belongs to, used to keep its name unique
   *  @return a Term holding `value` from `n` cycles ago
   */
  smt::Term make_history_chain(const smt::Term & value,
                               uint32_t n,
                               const std::string & name_prefix);

  /** Like make_history_chain(), but for a Bool-sorted (not BV-sorted)
   *  `cond` -- some solvers (e.g. Bitwuzla) can't build a state
   *  variable or a zero constant of sort Bool, so this wraps `cond`
   *  into a BV1 flag before delaying it and unwraps back to Bool
   *  afterward, matching the pattern already used for the `|=>`/
   *  `|-> ##N` antecedent delay.
   *  @param cond the current-cycle Bool-sorted condition to delay
   *  @param n    the number of cycles of delay (n == 0 is a no-op,
   *         returning `cond` unchanged)
   *  @param name_prefix see make_history_chain()
   *  @return a Bool-sorted Term equal to `cond` from `n` cycles ago
   */
  smt::Term delay_bool(const smt::Term & cond,
                       uint32_t n,
                       const std::string & name_prefix);

  /** Lazily create the shared "first cycle" flag: a 1-bit state var
   *  that is 1 in the initial state and 0 forever after.  Returned as
   *  a Boolean term (`flag == 1`).  Used to gate each LTL property's
   *  time-0 obligation so it is only asserted at the start of the
   *  trace.
   *  @param name_prefix used to name the flag the first time this is
   *         called (i.e. for whichever module's assertion first needs
   *         it) -- ignored on every subsequent call, since the flag
   *         is created once and cached for the whole design.
   */
  smt::Term init_flag(const std::string & name_prefix);

  /** Returns a Boolean term true iff the current cycle is one of the
   *  first `k` cycles of the trace (cycle index 0..k-1).  Used to
   *  exempt a `##k seq |-> ...` property from being checked before
   *  cycle `k`, the earliest cycle any `##k` match can end at (a
   *  naive per-cycle check would otherwise wrongly evaluate the
   *  consequent -- e.g. a `$past` with no real history yet -- at
   *  cycles no valid match could ever reach).
   *  @param k number of leading cycles to flag (k >= 1)
   *  @param name_prefix see make_history_chain()
   */
  smt::Term before_cycle(uint32_t k, const std::string & name_prefix);

  /** Returns a Boolean term true iff `disable_cond` holds at the
   *  current cycle or at any of the `window` cycles before it, or a
   *  null Term if `disable_cond` is null (no `disable iff` in scope).
   *  Used to extend a `disable iff` condition sampled at an
   *  antecedent's trigger cycle across the cycles a delayed (`##k`)
   *  consequent shifts the check over, so the whole match window is
   *  exempted, not just its last cycle.
   *  @param disable_cond the `disable iff` condition (explicit on the
   *         current assert statement, or the enclosing module's
   *         `default disable iff`) as a Boolean SMT term, or null
   *  @param window number of cycles before the current one to OR in
   *  @param name_prefix see make_history_chain()
   */
  smt::Term disable_window(const smt::Term & disable_cond,
                           uint32_t window,
                           const std::string & name_prefix);

  /** Tableau gadget builders.  Each takes the already-compiled
   *  `sat(...)` Boolean term(s) of the operand(s) and returns the
   *  `sat(...)` term of the composite temporal formula, installing the
   *  necessary "promise" latch and consistency constraint into the
   *  transition system.  The F / U builders also push their
   *  eventuality-discharge term onto `justice`.
   *  @param name_prefix see make_history_chain()
   */
  smt::Term make_X(const smt::Term & phi, const std::string & name_prefix);
  smt::Term make_G(const smt::Term & phi, const std::string & name_prefix);
  smt::Term make_F(const smt::Term & phi,
                   smt::TermVec & justice,
                   const std::string & name_prefix);
  smt::Term make_R(const smt::Term & a,
                   const smt::Term & b,
                   const std::string & name_prefix);
  smt::Term make_U(const smt::Term & a,
                   const smt::Term & b,
                   smt::TermVec & justice,
                   const std::string & name_prefix);

  /** Returns a fresh, monotonically increasing ID, drawn from the same
   *  counter as every latch this class itself mints. A caller minting
   *  its own hidden latch closely tied to one of this class's gadgets
   *  (e.g. a per-property activation latch alongside an LTL tableau
   *  obligation) uses this so its latch's name doesn't collide with
   *  (or reuse a number already used by) one of this class's own.
   */
  uint32_t next_id() { return latch_counter_++; }

 private:
  /** Build a hierarchical name for a latch: `name_prefix` (a module's
   *  hierarchical name prefix, or empty for the top level) followed by
   *  `name`.
   */
  std::string make_name(const std::string & name_prefix,
                        const std::string & name) const;

  FunctionalTransitionSystem & fts_;
  const smt::SmtSolver & solver_;

  // Monotonic counter used to mint unique state-var names for every
  // hidden latch this class introduces: history chains, the LTL
  // tableau's promise latches, and the before_cycle()/disable_window()
  // cycle-counting latches.
  uint32_t latch_counter_ = 0;

  // Cached Boolean term (`flag == 1`) for the shared LTL "first
  // cycle" state var, created on demand by init_flag().  Null until
  // the first LTL property is encoded.
  smt::Term ltl_init_flag_;
};

}  // namespace pono
