/*!
 * \file assertion_walker.h
 * \brief SVA/LTL AST-walking and dispatch: $past, sequences, and the
 *        LTL tableau's leaf/operator dispatch.
 * \author Áron Ricardo Perez-Lopez
 * \date 2026
 * \copyright See the LICENSE file in the top-level source directory.
 *
 * AssertionWalker owns everything to do with compiling an SVA assertion
 * statement into a safety property, a standing constraint, or an LTL
 * justice set: the concurrent/immediate assertion statement dispatch
 * (formerly two StatementKind cases in statement.cpp), the current-cycle-
 * Boolean fast path and general LTL tableau dispatch (formerly sva.cpp),
 * and the extracted safety properties (propvec()) and liveness justice
 * sets (ltl_justice()) themselves.
 *
 * Depends only on ExprEncoder (to convert the plain Boolean expressions at
 * the leaves of a property/sequence tree) and Tableau (for the tableau's
 * latch-building gadgets) -- both already independent of the rest of the
 * encoder, so this class is too: it holds no reference back to
 * SystemVerilogEncoder. `disable iff`'s enclosing-scope lookup
 * (Compilation::getDefaultDisable) is the caller's job -- process_
 * concurrent_assertion() takes the already-resolved default-disable
 * expression (or null) as a plain parameter instead of needing a
 * Compilation/Scope of its own.
 */
#pragma once

#include <string>
#include <vector>

#include "core/fts.h"
#include "smt-switch/smt.h"

namespace slang::ast {
class AssertionExpr;
class ConcurrentAssertionStatement;
class Expression;
class ImmediateAssertionStatement;
class Statement;
class Symbol;
class TimingControl;
}  // namespace slang::ast

namespace pono {

class ExprEncoder;
class Tableau;

class AssertionWalker
{
 public:
  AssertionWalker(ExprEncoder & expr_encoder,
                  Tableau & tableau,
                  const smt::SmtSolver & solver,
                  FunctionalTransitionSystem & fts);

  /** Process a concurrent assertion statement (`assert`/`assume`/
   *  `restrict`/`cover property (...)`, or `expect (...)`, which is
   *  logged and skipped as a simulation-only construct). Extends
   *  propvec()/ltl_justice() or adds a standing constraint to `fts`, as
   *  appropriate for the assertion kind and whether the property
   *  reduces to a current-cycle Boolean or needs the general LTL
   *  tableau.
   *  @param ca the concurrent assertion statement
   *  @param stmt the same statement, for its source label (used in log
   *         messages and in check_clock()'s exception message)
   *  @param prefix the caller's current hierarchical name prefix
   *  @param default_disable_expr the enclosing module's `default
   *         disable iff` condition (Compilation::getDefaultDisable),
   *         used only if `ca` has no explicit `disable iff` of its own;
   *         null if none applies
   */
  void process_concurrent_assertion(
      const slang::ast::ConcurrentAssertionStatement & ca,
      const slang::ast::Statement & stmt,
      const std::string & prefix,
      const slang::ast::Expression * default_disable_expr);

  /** Process a procedural immediate assertion (`assert`/`assume`/
   *  `restrict`/`cover (expr);`), guarded by the accumulated path
   *  `condition` (e.g. an enclosing `if`) rather than treated as
   *  always-active.
   *  @param ia the immediate assertion statement
   *  @param condition accumulated path condition (for if/case nesting)
   *  @param prefix the caller's current hierarchical name prefix
   */
  void process_immediate_assertion(
      const slang::ast::ImmediateAssertionStatement & ia,
      const smt::Term & condition,
      const std::string & prefix);

  /** @return the vector of safety properties (negated assertions)
   *  found so far.
   */
  smt::TermVec & propvec() { return propvec_; }

  /** @return the per-property generalized-Büchi justice sets found so
   *  far -- see SystemVerilogEncoder::Result::ltl_justice for what
   *  these mean.
   */
  std::vector<smt::TermVec> & ltl_justice() { return ltl_justice_; }

 private:
  /** Compile an SVA AssertionExpr into a Boolean SMT term that holds
   *  iff the assertion passes at the current cycle.  Returns a null
   *  Term when the expression uses an unsupported operator
   *  (e.g. liveness, sequence delays inside arbitrary positions,
   *  etc.); the caller can then skip that assertion.
   *  Non-overlapping implication (`|=>`) and the
   *  `|-> ##N` pattern introduce hidden latch state vars so the
   *  "P held N cycles ago" predicate is current-state-only.
   *  @param ae the assertion expression to compile
   *  @param prefix the current hierarchical name prefix
   *  @return the boolean term, or a null Term when unsupported
   */
  smt::Term assertion_expr_to_bool(const slang::ast::AssertionExpr & ae,
                                   const std::string & prefix);

  /** General bounded sequence matching: given a sequence expression
   *  (`Simple`/`SequenceWithMatch` with a consecutive `[m:n]`
   *  repetition, `SequenceConcat` with per-element `[m:n]` delay
   *  ranges, `FirstMatch`, `Clocking` -- checked via check_clock() and
   *  then unwrapped -- or a `Binary` intersect/within/throughout
   *  composition of two such sequences), returns a vector indexed by
   *  relative offset `L` where entry `L` is a Term true iff the
   *  sequence
   *  completes a match at the *current* cycle, having started `L`
   *  cycles earlier. A null entry means that offset is structurally
   *  unreachable. Returns an empty vector for sequence shapes this
   *  primitive doesn't (yet) model -- the caller should treat that the
   *  same as an unsupported construct.
   *
   *  Scoped to statically-bounded sequences: an unbounded (`[*]`,
   *  `[+]`, `[*n:$]`) or nonconsecutive/goto repetition, or an
   *  unbounded inter-element delay (`##[m:$]`), throws a clear
   *  PonoException rather than silently mismodeling or dropping it --
   *  this is a permanent architectural boundary of the encoder's
   *  compile-time-bounded model, not a "not implemented yet" gap.
   *  @param seq the sequence expression to match
   *  @param prefix the current hierarchical name prefix
   *  @return offsets indexed by relative start-to-end span
   */
  smt::TermVec offsets_ending_now(const slang::ast::AssertionExpr & seq,
                                  const std::string & prefix);

  /** Convenience wrapper over offsets_ending_now(): ORs together every
   *  reachable offset, i.e. "does `seq` complete a match at the
   *  current cycle, regardless of how long it took". Returns a null
   *  Term if offsets_ending_now() returns no reachable offsets at all
   *  (an unsupported sequence shape) so callers can fall back to their
   *  existing unsupported-construct handling.
   */
  smt::Term match_exists(const slang::ast::AssertionExpr & seq,
                         const std::string & prefix);

  /** The Boolean condition of a sequence's own leading element --
   *  "has an attempt to match `seq` just begun" -- used by
   *  weak_seq_bool() to detect when an in-progress match attempt has
   *  definitely failed. Recurses through FirstMatch/Clocking (like
   *  offsets_ending_now()) and into a SequenceConcat's first element.
   *  Throws for any other sequence shape (its own leading repetition,
   *  a `SequenceWithMatch`, or a `Binary` intersect/within/throughout
   *  as the outermost sequence) rather than guessing.
   */
  smt::Term leading_condition(const slang::ast::AssertionExpr & seq,
                              const std::string & prefix);

  /** Builds the `weak(seq)` Boolean safety condition: `seq` carries no
   *  obligation to ever match, but if an attempt began exactly
   *  `S = offsets_ending_now(seq).size() - 1` cycles ago (`S` being
   *  the sequence's own maximum span -- the last possible chance for
   *  that attempt to complete) and no completion happened anywhere in
   *  the intervening window, that attempt has definitely failed.
   *  Checked at every cycle, this covers every possible attempt start
   *  point exactly once, `S` cycles after it began. Returns a null
   *  Term if `seq`'s shape isn't modeled by offsets_ending_now().
   */
  smt::Term weak_seq_bool(const slang::ast::AssertionExpr & seq,
                          const std::string & prefix);

  /** General symbolic-tableau translation of an SVA property into the
   *  Boolean SMT term `sat(psi)` that holds at a cycle iff the
   *  (possibly negated) property `psi` holds starting from that cycle,
   *  where `psi` is `ae` when `neg` is false and `!ae` when `neg` is
   *  true.  Negation is pushed through the operators on the fly (so
   *  the gadgets built always correspond to the operators of `psi` in
   *  negation-normal form) -- this keeps the eventuality-fairness
   *  conditions correct regardless of the surrounding polarity.
   *
   *  Each temporal operator instantiates a one-step "promise" latch
   *  (see `tableau_`'s `make_X/G/F/R/U`) via `assign_next` plus a
   *  current-cycle consistency constraint, and every strong-
   *  eventuality operator (F / strong-until) appends its discharge
   *  condition to `justice`.
   *
   *  Returns a null Term when the property uses an operator the
   *  tableau does not model (sequence intersect/throughout/within/
   *  followed-by, etc.); the caller then skips the assertion.
   */
  smt::Term ltl_to_sat(const slang::ast::AssertionExpr & ae,
                       bool neg,
                       smt::TermVec & justice,
                       const std::string & prefix);

  /** Build the hierarchical name `prefix + "." + name` (or just `name`
   *  if `prefix` is empty) -- used for one hidden latch name and a few
   *  log messages; not worth a SymbolTable dependency for this trivial,
   *  stateless string concatenation (Tableau and SymbolTable each keep
   *  their own copy of this same one-liner for the same reason).
   */
  static std::string make_name(const std::string & prefix,
                               const std::string & name);

  /** Checks the clock named by a `Clocking` AssertionExpr node's
   *  `clocking` control against `design_clock_sym_`/`design_clock_edge_`,
   *  the (signal, edge) pair established by the first such clocking
   *  event seen anywhere in the design's properties. The first call
   *  overall just establishes that baseline; every later call throws a
   *  clear PonoException if it names a different signal or a different
   *  edge of the same signal -- this encoder has no clock-domain-
   *  crossing model (no clock dividers, no nondeterministic per-cycle
   *  choice of which clock toggles), so a genuinely multi-clock design
   *  is rejected outright rather than silently (or even just with a
   *  warning) collapsed onto one global cycle. See the "SVA design
   *  decisions" note at the top of assertion_walker.cpp. Also throws
   *  if `clocking` isn't a single edge-sensitive signal (`@*`, an
   *  event list, `repeat`, ...), since this check can't be soundly
   *  skipped for a shape it can't identify a clock from.
   *  @param clocking the clocking control to inspect
   */
  void check_clock(const slang::ast::TimingControl & clocking);

  ExprEncoder & expr_encoder_;
  Tableau & tableau_;
  const smt::SmtSolver & solver_;
  FunctionalTransitionSystem & fts_;

  // Safety properties extracted from SVA assert statements.
  smt::TermVec propvec_;

  // Per-property generalized-Büchi justice sets extracted from
  // temporal (LTL) assertions that are not pure safety.  Each entry
  // is the justice set { j_0, ..., j_k } of one assertion: a
  // counterexample is a lasso along which every j_i holds infinitely
  // often.  See SystemVerilogEncoder::Result::ltl_justice.
  std::vector<smt::TermVec> ltl_justice_;

  // The `disable iff` condition (explicit on the current assert
  // statement, or the enclosing module's `default disable iff`) as a
  // Boolean SMT term, or null if none applies.  Set just before
  // compiling one assertion's property expression and read by
  // assertion_expr_to_bool()/ltl_to_sat() via tableau_.disable_window(),
  // so it need not be threaded through every recursive call in between.
  smt::Term current_disable_cond_;

  // The (signal, edge) pair established by the first clocking event
  // seen anywhere in the design's properties (see check_clock()) --
  // persists for the lifetime of this AssertionWalker (one whole
  // design), not just one property, since the design has exactly one
  // clock or none at all. design_clock_edge_ is a slang::ast::EdgeKind
  // value stored as a plain int so this header doesn't need that
  // enum's full definition; only meaningful once design_clock_sym_ is
  // non-null.
  const slang::ast::Symbol * design_clock_sym_ = nullptr;
  int design_clock_edge_ = 0;

  // The current property's source label, for check_clock()'s exception
  // message; set at the start of each process_concurrent_assertion()
  // call. Immediate assertions have no AssertionExpr/Clocking tree to
  // walk, so process_immediate_assertion() never touches this.
  std::string current_assertion_label_;
};

}  // namespace pono
