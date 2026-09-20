/*!
 * \file assertion_walker.cpp
 * \brief SVA/LTL AST-walking and dispatch: $past, sequences, and the
 *        LTL tableau's leaf/operator dispatch.
 * \author Áron Ricardo Perez-Lopez
 * \date 2026
 * \copyright See the LICENSE file in the top-level source directory.
 *
 * See assertion_walker.h for what this class covers and why it holds no
 * reference back to SystemVerilogEncoder.
 *
 * SVA design decisions
 * ---------------------
 * A few SVA constructs have no single "obviously correct" encoding given
 * this encoder's model (a single global clock, and a propvec()-of-safety-
 * properties / ltl_justice()-of-liveness-obligations interface with no
 * notion of coverage or multiple clock domains). The choices made here
 * are deliberate and documented so future changes don't accidentally drift
 * from them:
 *
 *   - `cover property (P)` / immediate `cover (P)`: modeled via
 *     reachability duality -- checked exactly like `assert property (!P)`
 *     (or `assert (!P)`), so a "violation" of that surrogate assertion is
 *     precisely "P was reached". This is the standard way to expose a
 *     coverage goal through a safety-property-only interface. Temporal/
 *     sequence-shaped cover goals are out of scope (negating a liveness
 *     obligation isn't a reachability check) and throw a clear error.
 *
 *   - Multiclock properties (a property whose sequence mentions more than
 *     one clock, or more than one edge of the same clock): this encoder
 *     has no clock-domain-crossing model at all -- no clock dividers, no
 *     nondeterministic per-cycle choice of which clock toggles -- every
 *     design in this frontend's test suite already implicitly assumes one
 *     global clock advances the whole design by one cycle per sample.
 *     Rather than silently (or even just with a warning) collapse a
 *     second clock onto that same global cycle, check_clock() rejects it
 *     outright: the first `@(edge clk)` seen anywhere in the design's
 *     properties establishes the one clock this design is allowed to
 *     have, and any later property clocked on a different signal or a
 *     different edge throws a clear PonoException instead.
 */
#include "frontends/systemverilog/assertion_walker.h"

#include <cstdint>
#include <string>

#include "frontends/systemverilog/ast_helpers.h"
#include "frontends/systemverilog/bit_utils.h"
#include "frontends/systemverilog/expr_encoder.h"
#include "frontends/systemverilog/tableau.h"
#include "slang/ast/Expression.h"
#include "slang/ast/Symbol.h"
#include "slang/ast/TimingControl.h"
#include "slang/ast/expressions/AssertionExpr.h"
#include "slang/ast/expressions/MiscExpressions.h"
#include "slang/ast/statements/MiscStatements.h"
#include "slang/ast/types/Type.h"
#include "slang/syntax/AllSyntax.h"
#include "smt-switch/smt.h"
#include "utils/exceptions.h"
#include "utils/logger.h"

using namespace smt;
using namespace std;

namespace pono {

namespace {

// Detect a SequenceConcat that we can interpret as a constant
// k-cycle delay applied to a single inner assertion expression
// (`##k Q`).  Returns (k, Q*) on success, std::nullopt otherwise.
std::optional<std::pair<uint32_t, const slang::ast::AssertionExpr *>>
match_const_delay_seq(const slang::ast::AssertionExpr & ae)
{
  using namespace slang::ast;
  if (ae.kind != AssertionExprKind::SequenceConcat) return std::nullopt;
  auto & sc = ae.as<SequenceConcatExpr>();
  if (sc.elements.size() != 1) return std::nullopt;
  auto & e = sc.elements[0];
  if (!e.delay.max || *e.delay.max != e.delay.min) return std::nullopt;
  return std::make_pair(e.delay.min, &*e.sequence);
}

// The same shape, but for a delay with no upper bound (`##[m:$] Q`).
// Returns (m, Q*): the inner expression must hold at some cycle m or
// more from here, which is an eventuality and so belongs to the
// tableau rather than to the bounded sequence matcher.
std::optional<std::pair<uint32_t, const slang::ast::AssertionExpr *>>
match_unbounded_delay_seq(const slang::ast::AssertionExpr & ae)
{
  using namespace slang::ast;
  if (ae.kind != AssertionExprKind::SequenceConcat) return std::nullopt;
  auto & sc = ae.as<SequenceConcatExpr>();
  if (sc.elements.size() != 1) return std::nullopt;
  auto & e = sc.elements[0];
  if (e.delay.max) return std::nullopt;
  return std::make_pair(e.delay.min, &*e.sequence);
}

// A named `sequence`/`property` declaration referenced by name (e.g.
// `assert property (p_check);`) binds as a SimpleAssertionExpr wrapping
// an AssertionInstanceExpression -- not a plain boolean Expression --
// so routing it through expr_to_term() throws "unsupported expression
// kind". Slang has already expanded the referenced item's body with
// its actual arguments substituted -- a reference to a formal is
// itself an expanded instance -- so `body` is exactly the
// AssertionExpr this encoder should recurse into, arguments or not.
// Returns nullptr if `e` isn't such a reference (the caller should
// fall back to its normal expr_to_term() path).
//
// A local variable or a recursive instantiation still throws: each
// needs a binding environment of its own, which pre-expansion does
// not supply, so there is nothing correct to recurse into.
const slang::ast::AssertionExpr * resolve_named_assertion_ref(
    const slang::ast::Expression & e)
{
  using namespace slang::ast;
  if (e.kind != ExpressionKind::AssertionInstance) return nullptr;
  auto & aie = e.as<AssertionInstanceExpression>();
  if (aie.isRecursiveProperty || !aie.localVars.empty()) {
    throw PonoException(
        "SystemVerilogEncoder: named sequence/property references with "
        "local variables or recursion are not supported");
  }
  return &aie.body;
}

// Returns the source label of a concurrent assertion statement (e.g. the
// `p1` in `p1: assert property (...)`), or "<unnamed>" if it has none.
std::string assertion_label(const slang::ast::Statement & stmt)
{
  if (auto * syntax = stmt.syntax) {
    if (auto * ca_syntax =
            syntax
                ->as_if<slang::syntax::ConcurrentAssertionStatementSyntax>()) {
      if (ca_syntax->label) {
        return std::string(ca_syntax->label->name.valueText());
      }
    }
  }
  return "<unnamed>";
}

}  // namespace

AssertionWalker::AssertionWalker(ExprEncoder & expr_encoder,
                                 Tableau & tableau,
                                 const smt::SmtSolver & solver,
                                 FunctionalTransitionSystem & fts)
    : expr_encoder_(expr_encoder), tableau_(tableau), solver_(solver), fts_(fts)
{
}

string AssertionWalker::make_name(const string & prefix, const string & name)
{
  if (prefix.empty()) return name;
  return prefix + "." + name;
}

void AssertionWalker::check_clock(const slang::ast::TimingControl & clocking)
{
  using namespace slang::ast;

  if (clocking.kind != TimingControlKind::SignalEvent) {
    throw PonoException(
        "SystemVerilogEncoder: property '" + current_assertion_label_
        + "' has a clocking event that isn't a single edge-sensitive "
          "signal (@*, an event list, repeat, ...) -- this encoder "
          "requires every property to be clocked on one edge of one "
          "signal");
  }
  auto & sec = clocking.as<SignalEventControl>();
  const Symbol * sym = find_lhs_base(sec.expr);
  if (!sym) {
    throw PonoException(
        "SystemVerilogEncoder: property '" + current_assertion_label_
        + "' has a clocking event whose clock signal could not be "
          "resolved");
  }

  if (!design_clock_sym_) {
    design_clock_sym_ = sym;
    design_clock_edge_ = static_cast<int>(sec.edge);
    return;
  }
  if (sym != design_clock_sym_
      || static_cast<int>(sec.edge) != design_clock_edge_) {
    throw PonoException(
        "SystemVerilogEncoder: property '" + current_assertion_label_
        + "' is clocked on " + string(toString(sec.edge)) + " of '"
        + string(sym->name) + "', but this design's clock was already "
          "established as "
        + string(toString(static_cast<EdgeKind>(design_clock_edge_))) + " of '"
        + string(design_clock_sym_->name)
        + "' by an earlier property -- multi-clock designs (distinct "
          "clocks, or distinct edges of the same clock) are not "
          "supported, since this encoder has no clock-domain-crossing "
          "model (no clock dividers, no nondeterministic per-cycle "
          "choice of which clock toggles)");
  }
}

// ============================================================================
// LTL tableau dispatch
// ============================================================================
//
// Properties that are not pure safety are translated with a standard
// symbolic LTL tableau (temporal testers), built by tableau_'s
// make_X/G/F/R/U (see tableau.h/.cpp).  ltl_to_sat() below pushes negation
// to the leaves on the fly, so the testers it asks tableau_ to build always
// match the negation-normal form of the (negated) property.

namespace {
// The largest total span (in cycles) offsets_ending_now() will build
// before giving up -- a defensive cap against a pathological/absurd
// bounded sequence, mirroring the MAX_ITERS-style caps used elsewhere
// in this encoder for other compile-time-unrolled constructs.
constexpr uint32_t MAX_SEQ_WINDOW = 256;
}  // namespace

smt::TermVec AssertionWalker::offsets_ending_now(
    const slang::ast::AssertionExpr & seq,
    const string & prefix,
    bool * admits_empty)
{
  using namespace slang::ast;

  // Reports an empty match to a caller that asked about one, and
  // refuses to hide it from a caller that did not. Composing an
  // empty match is the concatenation rule's job; every other caller
  // would have to drop the alternative to carry on.
  auto report_empty = [&](const char * what) {
    if (admits_empty) {
      *admits_empty = true;
      return;
    }
    throw PonoException(string("SystemVerilogEncoder: ") + what
                        + " can match emptily, and an empty match occupies "
                          "no cycles, so there is no cycle at which this "
                          "position could say it completed");
  };

  // A single Boolean expression, optionally with a consecutive
  // repetition (`expr[*n:m]`, `expr[+]`, `expr[*]`). Shared by both
  // SimpleAssertionExpr and SequenceWithMatchExpr, which each carry
  // their own std::optional<SequenceRepetition>.
  auto boolean_with_repetition =
      [&](const slang::ast::Expression & expr,
          const std::optional<SequenceRepetition> & repetition) -> TermVec {
    Term b = expr_encoder_.expr_to_bool(expr, prefix);
    if (!repetition) return { b };
    if (repetition->kind != SequenceRepetition::Consecutive) {
      // Counting occurrences that need not be adjacent spans no
      // finite window. Decline, and ltl_to_sat() reads it as the
      // eventuality it is -- but only where a match merely has to
      // exist ahead, which the implication case below enforces.
      return {};
    }
    uint32_t lo = repetition->range.min;
    if (lo == 0) {
      // The empty alternative, which has no slot below. Reporting it
      // separately leaves exactly `[*1:hi]` behind, so `[*0:$]` needs
      // nothing beyond what `[*1:$]` already does.
      report_empty("a consecutive repetition with a zero lower bound");
      lo = 1;
    }
    // A run of at least `lo` ending now ends with a run of exactly
    // `lo`, and where the run began changes nothing about whether it
    // ends here -- so for these offsets the unbounded form and
    // `[*lo]` say the same thing.
    uint32_t hi = repetition->range.max ? *repetition->range.max : lo;
    if (hi >= MAX_SEQ_WINDOW) {
      throw PonoException("SystemVerilogEncoder: sequence repetition exceeds "
                          + std::to_string(MAX_SEQ_WINDOW) + " cycles");
    }
    // expr[*lo:hi] matches (ending now, started L cycles ago) for
    // L in [lo-1, hi-1], requiring expr to hold at every one of the
    // L+1 cycles from L-cycles-ago through now.
    TermVec out(hi, Term());
    Term running = b;
    for (uint32_t count = 1; count <= hi; ++count) {
      if (count > 1) {
        running = solver_->make_term(
            And, running, tableau_.make_history_chain(b, count - 1, prefix));
      }
      if (count >= lo) out[count - 1] = running;
    }
    return out;
  };

  // ORs together every non-null entry of a TermVec (the "does it
  // complete here at all" merge, without recomputing offsets_ending_now
  // -- used by Within below).
  auto or_vec = [&](const TermVec & v) -> Term {
    Term result;
    for (auto & t : v) {
      if (!t) continue;
      result = result ? solver_->make_term(Or, result, t) : t;
    }
    return result;
  };

  // OR of `base` and its delayed copies over the last `k` cycles
  // (`base` itself, plus 1..k cycles ago) -- "did `base` hold at *some*
  // point in the last k+1 cycles". Used by Within.
  auto window_or = [&](const Term & base, uint32_t k) -> Term {
    Term result = base;
    for (uint32_t j = 1; j <= k; ++j) {
      result = solver_->make_term(
          Or, result, tableau_.make_history_chain(base, j, prefix));
    }
    return result;
  };

  // AND of `base` and its delayed copies over the last `k` cycles --
  // "did `base` hold at *every* point in the last k+1 cycles". Used by
  // Throughout.
  auto window_and = [&](const Term & base, uint32_t k) -> Term {
    Term result = base;
    for (uint32_t j = 1; j <= k; ++j) {
      result = solver_->make_term(
          And, result, tableau_.make_history_chain(base, j, prefix));
    }
    return result;
  };

  switch (seq.kind) {
    case AssertionExprKind::Simple: {
      auto & simple = seq.as<SimpleAssertionExpr>();
      return boolean_with_repetition(simple.expr, simple.repetition);
    }

    case AssertionExprKind::SequenceWithMatch: {
      auto & swm = seq.as<SequenceWithMatchExpr>();
      // A parenthesized sequence with its own repetition
      // (`(seq)[*n:m]`) would need to convolve `seq`'s own offset
      // vector with itself count times -- out of scope; only a plain
      // Boolean operand with repetition is handled today.
      if (swm.repetition) {
        if (swm.expr.kind != AssertionExprKind::Simple
            || swm.expr.as<SimpleAssertionExpr>().repetition) {
          throw PonoException(
              "SystemVerilogEncoder: repetition of a non-Boolean sequence "
              "is not supported");
        }
        return boolean_with_repetition(swm.expr.as<SimpleAssertionExpr>().expr,
                                       swm.repetition);
      }
      return offsets_ending_now(swm.expr, prefix, admits_empty);
    }

    case AssertionExprKind::FirstMatch:
      // first_match(seq) only restricts *which* match is reported when
      // a sequence can match in more than one way -- it never changes
      // whether a match exists at all, which is all this encoder's
      // callers (an implication antecedent, an intersect/within/
      // throughout operand) ever ask offsets_ending_now() for.
      return offsets_ending_now(
          seq.as<FirstMatchAssertionExpr>().seq, prefix, admits_empty);

    case AssertionExprKind::Clocking: {
      // check_clock() throws if a nested clocking change inside a
      // sequence element names a different clock than the property's
      // first one -- see this file's multiclock design decision.
      auto & clk_expr = seq.as<ClockingAssertionExpr>();
      check_clock(clk_expr.clocking);
      return offsets_ending_now(clk_expr.expr, prefix, admits_empty);
    }

    case AssertionExprKind::SequenceConcat: {
      auto & sc = seq.as<SequenceConcatExpr>();
      TermVec acc;
      // Whether everything accumulated so far also matches emptily.
      // An empty match imposes no condition, so there is nothing to
      // store beyond the fact that it is available.
      bool acc_empty = false;

      // ORs `t` into `v` at index `idx`, growing `v` to fit. The
      // empty-match cases below land at smaller indices than the
      // plain ones, so the reachable width is easier to discover this
      // way than to compute up front.
      auto emit = [&](TermVec & v, size_t idx, const Term & t) {
        if (idx >= MAX_SEQ_WINDOW) {
          throw PonoException("SystemVerilogEncoder: sequence window exceeds "
                              + std::to_string(MAX_SEQ_WINDOW) + " cycles");
        }
        if (idx >= v.size()) v.resize(idx + 1, Term());
        v[idx] = v[idx] ? solver_->make_term(Or, v[idx], t) : t;
      };

      for (size_t i = 0; i < sc.elements.size(); ++i) {
        auto & elem = sc.elements[i];
        if (!elem.delay.max) {
          // No finite window spans this, so there are no offsets to
          // report; ltl_to_sat() takes it from here.
          return {};
        }
        uint32_t dmin = elem.delay.min;
        uint32_t dmax = *elem.delay.max;
        bool elem_empty = false;
        TermVec elem_offsets =
            offsets_ending_now(*elem.sequence, prefix, &elem_empty);
        // No offsets *and* no empty match is an unmodeled shape; no
        // offsets with one is `b[*0]`, which matches only emptily.
        if (elem_offsets.empty() && !elem_empty) return {};

        if (i == 0) {
          // The delay before the very first element just relabels how
          // far back "the sequence's start" is, with no extra
          // condition to AND in.
          for (size_t l = 0; l < elem_offsets.size(); ++l) {
            if (!elem_offsets[l]) continue;
            for (uint32_t d = dmin; d <= dmax; ++d) {
              emit(acc, l + d, elem_offsets[l]);
            }
          }
          if (elem_empty) {
            // A leading delay is genuinely `d` cycles of run-up, so
            // an empty first element leaves just those cycles -- and
            // with no delay either, nothing at all.
            for (uint32_t d = dmin; d <= dmax; ++d) {
              if (d == 0) {
                acc_empty = true;
              } else {
                emit(acc, d - 1, solver_->make_term(true));
              }
            }
          }
          continue;
        }

        TermVec new_acc;
        bool new_acc_empty = false;
        for (uint32_t d = dmin; d <= dmax; ++d) {
          // An inter-element delay counts from the prefix's last
          // cycle, so with no such cycle the LRM spends one less of
          // it: `(empty ##n s)` is `##(n-1) s`, and `(s ##n empty)`
          // is `s ##(n-1) 1`. Both collapse to `d - 1` run-up
          // cycles, or none at all when the delay is already zero.
          uint32_t dd = d == 0 ? 0 : d - 1;

          for (size_t lp = 0; lp < acc.size(); ++lp) {
            if (!acc[lp]) continue;
            for (size_t le = 0; le < elem_offsets.size(); ++le) {
              if (!elem_offsets[le]) continue;
              // Bring the (already-anchored-at-"now") prefix condition
              // back by (d + le) cycles so it aligns with this
              // element's own occurrence, then AND with this
              // element's own (unshifted) completion condition.
              Term shifted_prefix =
                  (d + le == 0)
                      ? acc[lp]
                      : tableau_.make_history_chain(acc[lp], d + le, prefix);
              emit(new_acc,
                   lp + d + le,
                   solver_->make_term(And, shifted_prefix, elem_offsets[le]));
            }
            if (elem_empty) {
              // The match ends where the prefix did, plus whatever
              // run-up cycles the delay still spends after it.
              Term shifted_prefix =
                  dd == 0 ? acc[lp]
                          : tableau_.make_history_chain(acc[lp], dd, prefix);
              emit(new_acc, lp + dd, shifted_prefix);
            }
          }

          if (acc_empty) {
            for (size_t le = 0; le < elem_offsets.size(); ++le) {
              if (!elem_offsets[le]) continue;
              // Nothing precedes this element, so it carries the
              // match on its own.
              emit(new_acc, dd + le, elem_offsets[le]);
            }
            if (elem_empty) {
              // Both sides empty: only the run-up cycles remain.
              if (dd == 0) {
                new_acc_empty = true;
              } else {
                emit(new_acc, dd - 1, solver_->make_term(true));
              }
            }
          }
        }
        acc = std::move(new_acc);
        acc_empty = new_acc_empty;
      }

      if (acc_empty) {
        report_empty("this sequence concatenation");
      }
      return acc;
    }

    case AssertionExprKind::Binary: {
      auto & b = seq.as<BinaryAssertionExpr>();
      switch (b.op) {
        case BinaryAssertionOperator::Intersect: {
          // s1 intersect s2 matches iff both match with the *same*
          // span -- AND the two offset vectors entry-by-entry.
          TermVec v1 = offsets_ending_now(b.left, prefix);
          TermVec v2 = offsets_ending_now(b.right, prefix);
          if (v1.empty() || v2.empty()) return {};
          size_t n = std::min(v1.size(), v2.size());
          TermVec out(n, Term());
          for (size_t k = 0; k < n; ++k) {
            if (v1[k] && v2[k]) out[k] = solver_->make_term(And, v1[k], v2[k]);
          }
          return out;
        }

        case BinaryAssertionOperator::Within: {
          // s1 within s2 matches over the same span as a match of s2,
          // provided s1 matches ending somewhere inside that span: for
          // each of s2's own completion offsets k, "s1 matched
          // somewhere in the last k+1 cycles" is exactly window_or()
          // over s1's merged "matches here" term.
          TermVec v1 = offsets_ending_now(b.left, prefix);
          TermVec v2 = offsets_ending_now(b.right, prefix);
          if (v1.empty() || v2.empty()) return {};
          Term s1_matches = or_vec(v1);
          if (!s1_matches) return {};
          TermVec out(v2.size(), Term());
          for (size_t k = 0; k < v2.size(); ++k) {
            if (!v2[k]) continue;
            out[k] = solver_->make_term(
                And, v2[k], window_or(s1_matches, (uint32_t)k));
          }
          return out;
        }

        case BinaryAssertionOperator::Throughout: {
          // expr throughout seq: the plain boolean expr must hold at
          // every cycle spanned by seq's match -- for each of seq's
          // own completion offsets k, that span is exactly the last
          // k+1 cycles, checked via window_and().
          Term expr_bool = assertion_expr_to_bool(b.left, prefix);
          if (!expr_bool) return {};
          TermVec v2 = offsets_ending_now(b.right, prefix);
          if (v2.empty()) return {};
          TermVec out(v2.size(), Term());
          for (size_t k = 0; k < v2.size(); ++k) {
            if (!v2[k]) continue;
            out[k] = solver_->make_term(
                And, v2[k], window_and(expr_bool, (uint32_t)k));
          }
          return out;
        }

        default:
          // And/Or/Iff/Implies/Until*/FollowedBy/etc. as a *sequence*
          // operand aren't sequence-composition operators in the SVA
          // sense (they combine boolean/property values, not match
          // spans) -- not something offsets_ending_now() is ever
          // asked for by this encoder's callers today.
          return {};
      }
    }

    default:
      // Unsupported sequence shape (a nested StrongWeak/etc. operand
      // this primitive doesn't model yet) -- the caller falls back to
      // its existing unsupported-construct handling.
      return {};
  }
}

smt::Term AssertionWalker::match_exists(const slang::ast::AssertionExpr & seq,
                                        const string & prefix)
{
  TermVec offsets = offsets_ending_now(seq, prefix);
  Term result;
  for (auto & t : offsets) {
    if (!t) continue;
    result = result ? solver_->make_term(Or, result, t) : t;
  }
  return result;
}

smt::Term AssertionWalker::leading_condition(
    const slang::ast::AssertionExpr & seq, const string & prefix)
{
  using namespace slang::ast;
  switch (seq.kind) {
    case AssertionExprKind::Simple: {
      auto & simple = seq.as<SimpleAssertionExpr>();
      if (simple.repetition) {
        // A repetition that can match emptily lets the attempt begin
        // without this expression holding at all, so the expression
        // no longer marks the start; every other consecutive count
        // needs its first iteration right here.
        if (simple.repetition->kind != SequenceRepetition::Consecutive
            || simple.repetition->range.min == 0) {
          throw PonoException(
              "SystemVerilogEncoder: weak()/strong() of a sequence whose "
              "leading repetition can match emptily is not supported");
        }
      }
      return expr_encoder_.expr_to_bool(simple.expr, prefix);
    }
    case AssertionExprKind::FirstMatch:
      return leading_condition(seq.as<FirstMatchAssertionExpr>().seq, prefix);
    case AssertionExprKind::Clocking: {
      auto & clk_expr = seq.as<ClockingAssertionExpr>();
      check_clock(clk_expr.clocking);
      return leading_condition(clk_expr.expr, prefix);
    }
    case AssertionExprKind::SequenceConcat:
      return leading_condition(
          *seq.as<SequenceConcatExpr>().elements[0].sequence, prefix);
    case AssertionExprKind::Binary: {
      // Which operand's start marks the whole sequence's start, read
      // off the same spans offsets_ending_now() builds for these.
      auto & b = seq.as<BinaryAssertionExpr>();
      switch (b.op) {
        case BinaryAssertionOperator::Intersect:
        case BinaryAssertionOperator::And:
          // Both operands start together (intersect also ends
          // together, which does not change where it begins).
          return solver_->make_term(And,
                                    leading_condition(b.left, prefix),
                                    leading_condition(b.right, prefix));
        case BinaryAssertionOperator::Or:
          return solver_->make_term(Or,
                                    leading_condition(b.left, prefix),
                                    leading_condition(b.right, prefix));
        case BinaryAssertionOperator::Within:
          // The span is the right operand's; the left may match
          // anywhere inside it, including later.
          return leading_condition(b.right, prefix);
        case BinaryAssertionOperator::Throughout: {
          // The span is the right operand's, and the left is a plain
          // expression that must hold at every cycle of it, so it
          // holds at the first one too.
          Term expr_bool = assertion_expr_to_bool(b.left, prefix);
          if (!expr_bool) {
            throw PonoException(
                "SystemVerilogEncoder: the left operand of `throughout` must "
                "be a plain expression");
          }
          return solver_->make_term(
              And, expr_bool, leading_condition(b.right, prefix));
        }
        default:
          throw PonoException(
              "SystemVerilogEncoder: weak()/strong() of this sequence "
              "operator is not supported");
      }
    }
    default:
      throw PonoException(
          "SystemVerilogEncoder: weak()/strong() of this sequence shape is "
          "not supported");
  }
}

smt::Term AssertionWalker::weak_seq_bool(const slang::ast::AssertionExpr & seq,
                                         const string & prefix)
{
  TermVec offsets = offsets_ending_now(seq, prefix);
  if (offsets.empty()) return Term();
  Term me;
  for (auto & t : offsets) {
    if (!t) continue;
    me = me ? solver_->make_term(Or, me, t) : t;
  }
  if (!me) return Term();

  // S = the sequence's own maximum span: the last possible cycle an
  // attempt that started here could still complete by.
  uint32_t s = static_cast<uint32_t>(offsets.size()) - 1;
  Term started_s_ago =
      tableau_.make_history_chain(leading_condition(seq, prefix), s, prefix);
  Term completed_in_window = me;
  for (uint32_t j = 1; j <= s; ++j) {
    completed_in_window = solver_->make_term(
        Or, completed_in_window, tableau_.make_history_chain(me, j, prefix));
  }
  // Violated iff an attempt began exactly S cycles ago and no
  // completion happened anywhere from then through now; weak(seq) is
  // the negation -- no obligation to ever attempt, but an attempt that
  // did begin must not be a definite, provable failure.
  Term violated = solver_->make_term(
      And, started_s_ago, solver_->make_term(Not, completed_in_window));
  return solver_->make_term(Not, violated);
}

namespace {

// Whether `ae` counts occurrences that need not be adjacent
// anywhere inside it. Such a count is an eventuality, which is the
// right reading for a consequent and the wrong one for an
// antecedent, so the implication case checks before recursing.
bool has_goto_repetition(const slang::ast::AssertionExpr & ae)
{
  using namespace slang::ast;
  switch (ae.kind) {
    case AssertionExprKind::Simple: {
      auto & simple = ae.as<SimpleAssertionExpr>();
      return simple.repetition
             && simple.repetition->kind != SequenceRepetition::Consecutive;
    }
    case AssertionExprKind::SequenceWithMatch: {
      auto & swm = ae.as<SequenceWithMatchExpr>();
      if (swm.repetition
          && swm.repetition->kind != SequenceRepetition::Consecutive) {
        return true;
      }
      return has_goto_repetition(swm.expr);
    }
    case AssertionExprKind::FirstMatch:
      return has_goto_repetition(ae.as<FirstMatchAssertionExpr>().seq);
    case AssertionExprKind::Clocking:
      return has_goto_repetition(ae.as<ClockingAssertionExpr>().expr);
    case AssertionExprKind::SequenceConcat: {
      for (auto & e : ae.as<SequenceConcatExpr>().elements) {
        if (has_goto_repetition(*e.sequence)) return true;
      }
      return false;
    }
    case AssertionExprKind::Binary: {
      auto & b = ae.as<BinaryAssertionExpr>();
      return has_goto_repetition(b.left) || has_goto_repetition(b.right);
    }
    default: return false;
  }
}

}  // namespace

smt::Term AssertionWalker::goto_repetition(
    const slang::ast::Expression & expr,
    const slang::ast::SequenceRepetition & rep,
    bool neg,
    smt::TermVec & justice,
    const string & prefix)
{
  using namespace slang::ast;
  if (rep.kind == SequenceRepetition::Consecutive) return Term();

  // `b[->n]` matches at the n-th occurrence of b; `b[=n]` may run on
  // past it, but its earliest match ends there too, and a consequent
  // only has to match somewhere -- so both come to the same thing
  // here.
  if (!rep.range.max || *rep.range.max != rep.range.min) {
    throw PonoException(
        "SystemVerilogEncoder: a goto or nonconsecutive repetition with a "
        "range of counts ([->m:n] / [=m:n]) is not supported");
  }
  uint32_t count = rep.range.min;
  if (count == 0) {
    throw PonoException(
        "SystemVerilogEncoder: a goto or nonconsecutive repetition of zero "
        "occurrences is not supported");
  }
  if (count > MAX_SEQ_WINDOW) {
    throw PonoException("SystemVerilogEncoder: repetition count exceeds "
                        + std::to_string(MAX_SEQ_WINDOW));
  }

  Term b = expr_encoder_.expr_to_bool(expr, prefix);
  if (!b) return Term();
  Term nb = solver_->make_term(Not, b);

  // Positively, the n-th occurrence is reached:
  //   P(1) = F b,  P(k) = F(b && X P(k-1))
  // Negated, that is its negation-normal form, which discharges no
  // eventuality and so emits no justice:
  //   N(1) = G !b, N(k) = G(!b || X N(k-1))
  Term acc =
      neg ? tableau_.make_G(nb, prefix) : tableau_.make_F(b, justice, prefix);
  for (uint32_t k = 2; k <= count; ++k) {
    Term next = tableau_.make_X(acc, prefix);
    acc = neg ? tableau_.make_G(solver_->make_term(Or, nb, next), prefix)
              : tableau_.make_F(
                    solver_->make_term(And, b, next), justice, prefix);
  }
  return acc;
}

smt::Term AssertionWalker::try_strong_sequence(
    const slang::ast::AssertionExpr & ae,
    bool neg,
    smt::TermVec & justice,
    const string & prefix)
{
  Term me = match_exists(ae, prefix);
  if (!me) return Term();
  if (in_weak_) {
    // Reached while unwrapping a weak() (see ltl_to_sat()'s
    // StrongWeak case). weak_seq_bool() is what models weak
    // correctly, and it only spans the shapes offsets_ending_now()
    // covers; getting here means it could not, so refuse rather than
    // attach the obligation weak explicitly withholds.
    throw PonoException(
        "SystemVerilogEncoder: weak() of this sequence shape is not "
        "supported -- it reduces to neither a bounded-span check nor a "
        "qualifier-free expression: "
        + current_assertion_label_);
  }
  // strong(seq): a genuine liveness obligation -- the sequence must
  // eventually complete a match.
  return neg ? tableau_.make_G(solver_->make_term(Not, me), prefix)
             : tableau_.make_F(me, justice, prefix);
}

smt::Term AssertionWalker::ltl_to_sat(const slang::ast::AssertionExpr & ae,
                                      bool neg,
                                      smt::TermVec & justice,
                                      const string & prefix)
{
  using namespace slang::ast;

  switch (ae.kind) {
    case AssertionExprKind::Clocking: {
      auto & clk_expr = ae.as<ClockingAssertionExpr>();
      check_clock(clk_expr.clocking);
      return ltl_to_sat(clk_expr.expr, neg, justice, prefix);
    }

    case AssertionExprKind::StrongWeak: {
      auto & sw = ae.as<StrongWeakAssertionExpr>();
      // The strong/weak qualifier only meaningfully differs for a
      // genuine bounded sequence -- must it eventually complete
      // (strong) or not (weak)? Any other shape (already a plain
      // Boolean/temporal expression) is unaffected by the qualifier
      // under this encoder's infinite-lasso semantics; just unwrap.
      // Unwrapping is right for a shape the qualifier cannot affect
      // (a plain Boolean or temporal expression), but any sequence
      // reached while unwrapping a *weak* one would pick up the
      // strong "must eventually complete" obligation from
      // try_strong_sequence() -- the opposite of what weak means, and
      // not necessarily at this node: in `weak(s1 and s2)` it is the
      // operands that acquire it. So the qualifier rides down the
      // recursion and try_strong_sequence() refuses under it.
      struct WeakScope
      {
        bool & flag;
        bool saved;
        ~WeakScope() { flag = saved; }
      } weak_scope{ in_weak_, in_weak_ };

      if (sw.strength == StrongWeakAssertionExpr::Strong) {
        // An explicit strong() inside a weak() is strong again.
        in_weak_ = false;
        if (Term strong = try_strong_sequence(sw.expr, neg, justice, prefix)) {
          return strong;
        }
      } else {
        in_weak_ = true;
      }
      return ltl_to_sat(sw.expr, neg, justice, prefix);
    }

    case AssertionExprKind::Simple: {
      auto & simple = ae.as<SimpleAssertionExpr>();
      if (auto * named = resolve_named_assertion_ref(simple.expr)) {
        return ltl_to_sat(*named, neg, justice, prefix);
      }
      if (simple.repetition) {
        if (Term counted = goto_repetition(
                simple.expr, *simple.repetition, neg, justice, prefix)) {
          return counted;
        }
        // See the matching check in assertion_expr_to_bool(): route
        // through the general bounded sequence matcher (which throws
        // for an unbounded repeat count) instead of silently ignoring
        // the repetition.
        Term me = match_exists(ae, prefix);
        if (!me) return Term();
        return neg ? solver_->make_term(Not, me) : me;
      }
      Term b = expr_encoder_.expr_to_bool(simple.expr, prefix);
      if (!b) return Term();
      return neg ? solver_->make_term(Not, b) : b;
    }

    case AssertionExprKind::SequenceConcat: {
      // A bare `##k Q` property has the same truth value as Q under
      // our infinite-time semantics (modulo a front shift).  Unwrap.
      if (auto m = match_const_delay_seq(ae)) {
        return ltl_to_sat(*m->second, neg, justice, prefix);
      }
      if (auto u = match_unbounded_delay_seq(ae)) {
        Term inner = ltl_to_sat(*u->second, neg, justice, prefix);
        if (!inner) return Term();
        inner = neg ? tableau_.make_G(inner, prefix)
                    : tableau_.make_F(inner, justice, prefix);
        for (uint32_t i = 0; i < u->first; ++i) {
          inner = tableau_.make_X(inner, prefix);
        }
        return inner;
      }
      // A genuine multi-element sequence used directly as a property
      // (`assert property (a ##1 b);`, as opposed to as the
      // antecedent of `|->`/`|=>`, which assertion_expr_to_bool()
      // handles via the bounded sequence matcher): per the LRM this
      // has implicit `strong` semantics -- the sequence must
      // eventually complete a match.
      if (Term strong = try_strong_sequence(ae, neg, justice, prefix)) {
        return strong;
      }
      throw PonoException(
          "SystemVerilogEncoder: property '" + current_assertion_label_
          + "' uses a sequence shape used directly as a property that "
            "is not supported");
    }

    case AssertionExprKind::Unary: {
      auto & u = ae.as<UnaryAssertionExpr>();
      if (u.op == UnaryAssertionOperator::Not) {
        return ltl_to_sat(u.expr, !neg, justice, prefix);
      }

      bool is_always = u.op == UnaryAssertionOperator::Always
                       || u.op == UnaryAssertionOperator::SAlways;
      bool is_next = u.op == UnaryAssertionOperator::NextTime
                     || u.op == UnaryAssertionOperator::SNextTime;
      // Every UnaryAssertionOperator value is handled here; this is
      // defensive against a future slang operator this encoder hasn't
      // been taught, not a currently-reachable case.
      if (!is_always && !is_next && u.op != UnaryAssertionOperator::Eventually
          && u.op != UnaryAssertionOperator::SEventually) {
        throw PonoException(
            "SystemVerilogEncoder: unsupported unary assertion operator "
            + string(toString(u.op)));
      }

      // Normalize the optional cycle window `[m:n]` so the unranged
      // forms fall out of the same construction: a bare `nexttime` is
      // `[1:1]`, and a bare `always`/`s_eventually` is `[0:$]`, which
      // reduces to the plain G/F tester below.  The weak/strong pairs
      // coincide here -- that distinction only bites at the end of a
      // truncated trace, and Pono reasons about infinite behaviours.
      uint32_t lo = u.range ? u.range->min : (is_next ? 1 : 0);
      std::optional<uint32_t> hi =
          u.range ? u.range->max
                  : (is_next ? std::optional<uint32_t>(1) : std::nullopt);
      if (hi && *hi >= MAX_SEQ_WINDOW) {
        throw PonoException(
            "SystemVerilogEncoder: property cycle range "
            "exceeds "
            + std::to_string(MAX_SEQ_WINDOW) + " cycles");
      }
      if (hi && *hi < lo) {
        throw PonoException(
            "SystemVerilogEncoder: property cycle range has a maximum "
            "below its minimum");
      }

      // ltl_to_sat() already returns sat(!p) when `neg`, and !X == X!,
      // so the negation rides along inside phi -- only the way the
      // per-cycle terms are combined flips.
      Term phi = ltl_to_sat(u.expr, neg, justice, prefix);
      if (!phi) return Term();

      if (!hi) {
        // `always [m:$]` / `s_eventually [m:$]` (the only two operators
        // slang lets go unbounded): shift the unbounded tester forward
        // by m.  m == 0 is the bare `always`/`s_eventually` form.
        Term inner = (is_always != neg) ? tableau_.make_G(phi, prefix)
                                        : tableau_.make_F(phi, justice, prefix);
        for (uint32_t i = 0; i < lo; ++i) {
          inner = tableau_.make_X(inner, prefix);
        }
        return inner;
      }

      // Bounded window: AND (always) or OR (eventually) of phi shifted
      // forward j cycles, for j in [lo, hi].  Built as one chain so the
      // shared X^lo..X^j prefix is not rebuilt per j -- that is `hi`
      // make_X() testers in total rather than O(hi^2).
      PrimOp combine = (is_always != neg) ? And : Or;
      Term shifted = phi;
      for (uint32_t i = 0; i < lo; ++i) {
        shifted = tableau_.make_X(shifted, prefix);
      }
      Term result = shifted;  // X^lo phi
      for (uint32_t j = lo + 1; j <= *hi; ++j) {
        shifted = tableau_.make_X(shifted, prefix);  // X^j phi
        result = solver_->make_term(combine, result, shifted);
      }
      return result;
    }

    case AssertionExprKind::Conditional: {
      // In-property `if (cond) p else q`: negation distributes into
      // whichever branch `cond` selects -- the branch structure itself
      // doesn't change polarity, so this is a plain ITE over the two
      // (already correctly negated) recursive results.
      auto & c = ae.as<ConditionalAssertionExpr>();
      Term cond_bool = expr_encoder_.expr_to_bool(c.condition, prefix);
      Term if_branch = ltl_to_sat(c.ifExpr, neg, justice, prefix);
      if (!if_branch) return Term();
      Term else_branch;
      if (c.elseExpr) {
        else_branch = ltl_to_sat(*c.elseExpr, neg, justice, prefix);
        if (!else_branch) return Term();
      } else {
        // No else branch: per the LRM, a false condition with nothing
        // to check is vacuously satisfied -- true when un-negated,
        // false negated.
        else_branch = solver_->make_term(!neg);
      }
      return solver_->make_term(Ite, cond_bool, if_branch, else_branch);
    }

    case AssertionExprKind::Case: {
      // In-property `case (sel) item0: p0; item1: p1; ... endcase`:
      // the same ITE idea as Conditional, generalized to N branches --
      // fold the item list right-to-left into a chain of ITEs seeded
      // by the default case's (already correctly negated) result.
      auto & c = ae.as<CaseAssertionExpr>();
      Term sel = expr_encoder_.expr_to_term(c.expr, prefix);
      Term result;
      if (c.defaultCase) {
        result = ltl_to_sat(*c.defaultCase, neg, justice, prefix);
        if (!result) return Term();
      } else {
        // No default and nothing matched: vacuously satisfied, same
        // as Conditional's missing else branch above.
        result = solver_->make_term(!neg);
      }
      uint64_t sel_w = sel->get_sort()->get_width();
      for (auto it = c.items.rbegin(); it != c.items.rend(); ++it) {
        Term branch = ltl_to_sat(*it->body, neg, justice, prefix);
        if (!branch) return Term();
        Term item_cond;
        for (auto * match_expr : it->expressions) {
          Term m = expr_encoder_.expr_to_term(*match_expr, prefix);
          m = resize_to(solver_, m, sel_w, match_expr->type->isSigned());
          Term eq = solver_->make_term(Equal, sel, m);
          item_cond = item_cond ? solver_->make_term(Or, item_cond, eq) : eq;
        }
        result = solver_->make_term(Ite, item_cond, branch, result);
      }
      return result;
    }

    case AssertionExprKind::Binary: {
      auto & b = ae.as<BinaryAssertionExpr>();
      switch (b.op) {
        case BinaryAssertionOperator::And:
        case BinaryAssertionOperator::Or: {
          Term l = ltl_to_sat(b.left, neg, justice, prefix);
          Term r = ltl_to_sat(b.right, neg, justice, prefix);
          if (!l || !r) return Term();
          bool is_and = (b.op == BinaryAssertionOperator::And);
          if (neg) is_and = !is_and;  // De Morgan
          return solver_->make_term(is_and ? And : Or, l, r);
        }

        case BinaryAssertionOperator::Iff: {
          // a iff b == (a implies b) && (b implies a); build both
          // polarities of each operand through the recursion
          // (unavoidable -- iff isn't monotone in either argument) so
          // a temporal operand gets its own properly
          // negation-normalized dual construction, the same as every
          // other case above, instead of only ever being visited
          // positively.
          Term na = ltl_to_sat(b.left, true, justice, prefix);
          Term a = ltl_to_sat(b.left, false, justice, prefix);
          Term nb = ltl_to_sat(b.right, true, justice, prefix);
          Term bb = ltl_to_sat(b.right, false, justice, prefix);
          if (!na || !a || !nb || !bb) return Term();
          Term a_implies_b = solver_->make_term(Or, na, bb);
          Term b_implies_a = solver_->make_term(Or, nb, a);
          if (!neg) return solver_->make_term(And, a_implies_b, b_implies_a);
          // !(a iff b) == (a && !b) || (b && !a)
          return solver_->make_term(Or,
                                    solver_->make_term(And, a, nb),
                                    solver_->make_term(And, bb, na));
        }

        case BinaryAssertionOperator::Implies: {
          // a implies b == !a || b ;  !(a implies b) == a && !b
          Term l = ltl_to_sat(b.left, !neg, justice, prefix);
          Term r = ltl_to_sat(b.right, neg, justice, prefix);
          if (!l || !r) return Term();
          return solver_->make_term(neg ? And : Or, l, r);
        }

        case BinaryAssertionOperator::OverlappedImplication:
        case BinaryAssertionOperator::NonOverlappedImplication: {
          // seq |-> prop / seq |=> prop with a Boolean antecedent:
          //   !a || X^delay b   (delay 0 overlapped, 1 non-overlapped,
          //                       plus any ##k on the consequent)
          // and its negation a && X^delay !b.
          uint32_t delay =
              (b.op == BinaryAssertionOperator::NonOverlappedImplication) ? 1
                                                                          : 0;
          const AssertionExpr * rhs = &b.right;
          // `##[m:$] Q` says Q holds at some cycle m or more away,
          // which is an eventuality: the delay fixes where the wait
          // starts, and F carries it the rest of the way.
          bool unbounded = false;
          if (auto m = match_const_delay_seq(b.right)) {
            delay += m->first;
            rhs = m->second;
          } else if (auto u = match_unbounded_delay_seq(b.right)) {
            delay += u->first;
            rhs = u->second;
            unbounded = true;
          }
          if (has_goto_repetition(b.left)) {
            // An antecedent has to say its match ends *now*, which
            // for a count of non-adjacent occurrences is a fact
            // about unbounded history, not an eventuality -- so the
            // reading ltl_to_sat() would give it is the wrong one.
            throw PonoException(
                "SystemVerilogEncoder: a goto or nonconsecutive repetition "
                "in an implication's antecedent is not supported");
          }
          Term l = ltl_to_sat(b.left, !neg, justice, prefix);
          Term r = ltl_to_sat(*rhs, neg, justice, prefix);
          if (!l || !r) return Term();
          if (unbounded) {
            // Negated, "never at or after m" is a safety obligation
            // and discharges nothing, so only the positive form
            // contributes a justice condition.
            r = neg ? tableau_.make_G(r, prefix)
                    : tableau_.make_F(r, justice, prefix);
          }
          for (uint32_t i = 0; i < delay; ++i) r = tableau_.make_X(r, prefix);
          return solver_->make_term(neg ? And : Or, l, r);
        }

        case BinaryAssertionOperator::Until:
        case BinaryAssertionOperator::SUntil:
        case BinaryAssertionOperator::UntilWith:
        case BinaryAssertionOperator::SUntilWith: {
          bool strong = (b.op == BinaryAssertionOperator::SUntil
                         || b.op == BinaryAssertionOperator::SUntilWith);
          bool with = (b.op == BinaryAssertionOperator::UntilWith
                       || b.op == BinaryAssertionOperator::SUntilWith);
          if (!neg) {
            Term l = ltl_to_sat(b.left, false, justice, prefix);
            Term r = ltl_to_sat(b.right, false, justice, prefix);
            if (!l || !r) return Term();
            // until_with: the terminating cycle must also satisfy a.
            Term term = with ? solver_->make_term(And, l, r) : r;
            if (strong) return tableau_.make_U(l, term, justice, prefix);
            // weak until: a W term == term R (a || term).
            return tableau_.make_R(
                term, solver_->make_term(Or, l, term), prefix);
          }
          // Negated, with operands already in negation-normal form
          // (nl = sat(!left), nr = sat(!right)):
          //   !(a U_strong b)      = !a R !b
          //   !(a W b)             = !b U (!a && !b)
          //   !(a U_strong (a&&b)) = !a R (!a || !b)
          //   !(a W (a&&b))        = (!a || !b) U !a
          Term nl = ltl_to_sat(b.left, true, justice, prefix);
          Term nr = ltl_to_sat(b.right, true, justice, prefix);
          if (!nl || !nr) return Term();
          if (!with) {
            if (strong) return tableau_.make_R(nl, nr, prefix);
            return tableau_.make_U(
                nr, solver_->make_term(And, nl, nr), justice, prefix);
          }
          Term nterm = solver_->make_term(Or, nl, nr);  // !(a && b)
          if (strong) return tableau_.make_R(nl, nterm, prefix);
          return tableau_.make_U(nterm, nl, justice, prefix);
        }

        case BinaryAssertionOperator::OverlappedFollowedBy:
        case BinaryAssertionOperator::NonOverlappedFollowedBy: {
          // s1 #-# p2 / s1 #=# p2 is a required (not merely
          // conditional, unlike |->) sequential composition: s1 must
          // match, and p2 must then hold starting at the match's own
          // end cycle (overlapped) or one cycle later (non-
          // overlapped). Since every property here is checked at
          // *every* cycle (assert property(p) ~ always p), reindex so
          // "now" is p2's own check point -- exactly the trick
          // OverlappedImplication/NonOverlappedImplication already
          // use above (delay the antecedent match *backward* to align
          // with the consequent's check point, rather than shifting
          // the consequent forward): match_exists(s1) already gives
          // "s1 matched, ending now" for any match length; delaying
          // that *backward* by `extra` cycles gives "s1 matched,
          // ending `extra` cycles ago" -- exactly the antecedent this
          // formula needs when p2 is checked now.
          Term match = match_exists(b.left, prefix);
          if (!match) return Term();
          uint32_t extra =
              (b.op == BinaryAssertionOperator::NonOverlappedFollowedBy) ? 1
                                                                         : 0;
          if (extra > 0)
            match = tableau_.make_history_chain(match, extra, prefix);
          Term p2 = ltl_to_sat(b.right, neg, justice, prefix);
          if (!p2) return Term();
          // sat(s1 #-# p2) = match AND p2; negated (De Morgan):
          // !match OR sat(!p2) -- p2 is already the correctly negated
          // term per this function's usual convention, so only
          // `match`'s own polarity and the outer combinator branch on
          // `neg`.
          Term match_term = neg ? solver_->make_term(Not, match) : match;
          return solver_->make_term(neg ? Or : And, match_term, p2);
        }

        default:
          // Intersect / Throughout / Within: sequence-composition
          // operators offsets_ending_now() already models when the
          // whole binary expression is treated as a sequence (see
          // SeqIntersect/SeqWithin/SeqThroughout) -- try that (the
          // same implicit-strong fallback SequenceConcat/FirstMatch/
          // SequenceWithMatch use below) before giving up.
          if (Term strong = try_strong_sequence(ae, neg, justice, prefix)) {
            return strong;
          }
          throw PonoException(
              "SystemVerilogEncoder: property '" + current_assertion_label_
              + "' uses '" + string(toString(b.op))
              + "' as a top-level connective, which is not supported");
      }
    }

    // FirstMatch, SequenceWithMatch: bounded sequence shapes
    // offsets_ending_now() already models -- try treating the whole
    // property as an implicitly-strong sequence match (the same
    // fallback the SequenceConcat/Binary cases above use) before
    // giving up. Abort (accept_on/reject_on/sync_accept_on/sync_
    // reject_on) and a nested DisableIff (one not stripped by the
    // top-level `disable iff` handling in
    // process_concurrent_assertion(), e.g. as one operand of a
    // Binary/Unary operator) aren't sequences at all, so this always
    // fails for them, falling through to the throw below.
    default:
      if (Term strong = try_strong_sequence(ae, neg, justice, prefix)) {
        return strong;
      }
      throw PonoException(
          "SystemVerilogEncoder: property '" + current_assertion_label_
          + "' uses an assertion expression shape (" + string(toString(ae.kind))
          + ") that is not supported");
  }
}

smt::Term AssertionWalker::reanchor(const Term & t,
                                    uint32_t from,
                                    uint32_t to,
                                    const string & prefix)
{
  return from == to ? t : tableau_.make_history_chain(t, to - from, prefix);
}

smt::Term AssertionWalker::assertion_expr_to_bool(
    const slang::ast::AssertionExpr & ae, const string & prefix)
{
  uint32_t span = 0;
  Term t = assertion_expr_to_bool(ae, prefix, span);
  // A re-anchored term is only meaningful once gated with
  // before_cycle(span).  Callers of this form have nowhere to put
  // that gate, so decline and let them fall back to the tableau; the
  // one caller that can gate uses the three-argument form.
  return span == 0 ? t : Term();
}

smt::Term AssertionWalker::assertion_expr_to_bool(
    const slang::ast::AssertionExpr & ae,
    const string & prefix,
    uint32_t & span)
{
  using namespace slang::ast;

  span = 0;

  switch (ae.kind) {
    case AssertionExprKind::Clocking: {
      // The clocking event has already been baked into our cycle
      // abstraction; check_clock() throws if it names a different
      // clock than the property's first one, then just recurse into
      // the underlying expression.
      auto & clk_expr = ae.as<ClockingAssertionExpr>();
      check_clock(clk_expr.clocking);
      return assertion_expr_to_bool(clk_expr.expr, prefix, span);
    }

    case AssertionExprKind::Simple: {
      auto & simple = ae.as<SimpleAssertionExpr>();
      if (auto * named = resolve_named_assertion_ref(simple.expr)) {
        return assertion_expr_to_bool(*named, prefix, span);
      }
      if (simple.repetition) {
        // `expr[*n:m]`/`expr[+]`/`expr[*]`: route through the general
        // bounded sequence matcher (which throws for an unbounded
        // repeat count) instead of silently ignoring the repetition
        // and returning a bare `bool(expr)`.
        return match_exists(ae, prefix);
      }
      return expr_encoder_.expr_to_bool(simple.expr, prefix);
    }

    case AssertionExprKind::SequenceConcat: {
      // A standalone `##k Q` as a property: in our infinite-time
      // safety encoding the constant front-shift doesn't change the
      // truth value (it just postpones when the first violation can
      // be reported), so unwrap the inner sequence.
      if (auto matched = match_const_delay_seq(ae)) {
        return assertion_expr_to_bool(*matched->second, prefix, span);
      }
      return Term();
    }

    case AssertionExprKind::StrongWeak: {
      auto & sw = ae.as<StrongWeakAssertionExpr>();
      if (sw.strength == StrongWeakAssertionExpr::Strong) {
        // strong(seq) is a genuine liveness obligation ("must
        // eventually complete"), not reducible to a current-cycle
        // Boolean -- handled by ltl_to_sat()'s StrongWeak case
        // instead, which builds the eventuality tableau. Returning
        // null here forces the ConcurrentAssertion handler to fall
        // through to that path rather than (incorrectly) treating a
        // strong sequence as if it were always-true right now.
        return Term();
      }
      // weak(seq): no obligation to ever match, but an attempt that
      // did begin must not be a definite, provable failure -- see
      // weak_seq_bool(). Any other shape (already a plain Boolean/
      // temporal expression) is unaffected by the qualifier; just
      // unwrap.
      if (Term w = weak_seq_bool(sw.expr, prefix)) return w;
      return assertion_expr_to_bool(sw.expr, prefix);
    }

    case AssertionExprKind::Unary: {
      auto & u = ae.as<UnaryAssertionExpr>();
      if (u.op == UnaryAssertionOperator::Not) {
        // Negation keeps the operand's anchor, so the span rides
        // through unchanged.
        Term inner = assertion_expr_to_bool(u.expr, prefix, span);
        if (!inner) return Term();
        return solver_->make_term(Not, inner);
      }

      bool is_always = u.op == UnaryAssertionOperator::Always
                       || u.op == UnaryAssertionOperator::SAlways;
      bool is_next = u.op == UnaryAssertionOperator::NextTime
                     || u.op == UnaryAssertionOperator::SNextTime;
      bool is_eventually = u.op == UnaryAssertionOperator::Eventually
                           || u.op == UnaryAssertionOperator::SEventually;
      if (!is_always && !is_next && !is_eventually) return Term();

      uint32_t lo = 0;
      uint32_t hi = 0;
      if (u.range) {
        // `[m:$]` reaches arbitrarily far forward, so there is no last
        // cycle to re-anchor to -- that one stays with the tableau.
        if (!u.range->max) return Term();
        lo = u.range->min;
        hi = *u.range->max;
      } else if (is_next) {
        // A bare `nexttime` is `[1:1]`.
        lo = hi = 1;
      } else if (is_always) {
        // Unranged `always P`: the whole-trace closure is already
        // implicit in the per-cycle property check, so P alone is it.
        return assertion_expr_to_bool(u.expr, prefix, span);
      } else {
        // Unranged `s_eventually`: unbounded, genuinely liveness.
        return Term();
      }
      if (hi < lo || hi >= MAX_SEQ_WINDOW) return Term();

      uint32_t inner_span = 0;
      Term inner = assertion_expr_to_bool(u.expr, prefix, inner_span);
      if (!inner) return Term();

      // Every cycle the window names is in the past once the check is
      // re-anchored to the window's last one, so the operand read `j`
      // cycles after the attempt started is read `hi - j` cycles ago.
      // `always` needs all of them, `eventually` any of them.
      span = hi + inner_span;
      Term result;
      for (uint32_t j = lo; j <= hi; ++j) {
        Term shifted = reanchor(inner, j + inner_span, span, prefix);
        result = result
                     ? solver_->make_term(is_always ? And : Or, result, shifted)
                     : shifted;
      }
      return result;
    }

    case AssertionExprKind::Conditional: {
      // Mirrors ltl_to_sat()'s Conditional case (always positive
      // polarity here; the caller negates the whole result if
      // needed), so a purely-Boolean if/else takes this fast,
      // current-cycle-safety path instead of unconditionally paying
      // for the full LTL tableau. Falls back (returns null) as soon
      // as either branch does.
      auto & c = ae.as<ConditionalAssertionExpr>();
      uint32_t if_span = 0;
      Term if_branch = assertion_expr_to_bool(c.ifExpr, prefix, if_span);
      if (!if_branch) return Term();
      uint32_t else_span = 0;
      Term else_branch;
      if (c.elseExpr) {
        else_branch = assertion_expr_to_bool(*c.elseExpr, prefix, else_span);
        if (!else_branch) return Term();
      } else {
        else_branch = solver_->make_term(true);
      }
      // Bring both branches, and the condition they select between,
      // onto whichever anchor is later.
      span = std::max(if_span, else_span);
      if_branch = reanchor(if_branch, if_span, span, prefix);
      else_branch = reanchor(else_branch, else_span, span, prefix);
      Term cond_bool = reanchor(
          expr_encoder_.expr_to_bool(c.condition, prefix), 0, span, prefix);
      return solver_->make_term(Ite, cond_bool, if_branch, else_branch);
    }

    case AssertionExprKind::Case: {
      // Mirrors ltl_to_sat()'s Case case -- see Conditional above.
      auto & c = ae.as<CaseAssertionExpr>();
      // Convert every branch before building anything: they all have
      // to agree on an anchor before the selector can choose between
      // them.
      uint32_t default_span = 0;
      Term default_branch;
      if (c.defaultCase) {
        default_branch =
            assertion_expr_to_bool(*c.defaultCase, prefix, default_span);
        if (!default_branch) return Term();
      } else {
        default_branch = solver_->make_term(true);
      }
      TermVec bodies;
      std::vector<uint32_t> body_spans;
      span = default_span;
      for (auto & item : c.items) {
        uint32_t body_span = 0;
        Term body = assertion_expr_to_bool(*item.body, prefix, body_span);
        if (!body) return Term();
        bodies.push_back(body);
        body_spans.push_back(body_span);
        span = std::max(span, body_span);
      }

      Term sel =
          reanchor(expr_encoder_.expr_to_term(c.expr, prefix), 0, span, prefix);
      uint64_t sel_w = sel->get_sort()->get_width();
      Term result = reanchor(default_branch, default_span, span, prefix);
      for (size_t i = c.items.size(); i-- > 0;) {
        Term branch = reanchor(bodies[i], body_spans[i], span, prefix);
        Term item_cond;
        for (auto * match_expr : c.items[i].expressions) {
          Term m = expr_encoder_.expr_to_term(*match_expr, prefix);
          m = resize_to(solver_, m, sel_w, match_expr->type->isSigned());
          Term eq =
              solver_->make_term(Equal, sel, reanchor(m, 0, span, prefix));
          item_cond = item_cond ? solver_->make_term(Or, item_cond, eq) : eq;
        }
        result = solver_->make_term(Ite, item_cond, branch, result);
      }
      return result;
    }

    case AssertionExprKind::Binary: {
      auto & b = ae.as<BinaryAssertionExpr>();
      bool is_impl =
          (b.op == BinaryAssertionOperator::OverlappedImplication
           || b.op == BinaryAssertionOperator::NonOverlappedImplication
           || b.op == BinaryAssertionOperator::Implies);

      if (is_impl) {
        // A `##k` prefix on the antecedent restricts which cycles a
        // match can even start from (the earliest anchor cycle is
        // k); gate the whole implication so it is vacuously true
        // before cycle k instead of dropping the delay, which would
        // otherwise evaluate the consequent (e.g. a `$past` with no
        // real history yet) at cycles no valid match could reach.
        uint32_t lhs_delay = 0;
        const AssertionExpr * lhs_inner = &b.left;
        if (auto lhs_matched = match_const_delay_seq(b.left)) {
          lhs_delay = lhs_matched->first;
          lhs_inner = lhs_matched->second;
        }
        // A plain-Boolean-reducible antecedent (the common case) is
        // handled directly; a multi-element/first-match/nested-clock
        // sequence antecedent (`a ##1 b |-> ...`,
        // `first_match(seq) |-> ...`) falls back to the general
        // bounded sequence matcher.
        uint32_t lhs_span = 0;
        Term lhs = assertion_expr_to_bool(*lhs_inner, prefix, lhs_span);
        if (!lhs) {
          lhs = match_exists(*lhs_inner, prefix);
          lhs_span = 0;
        }
        if (!lhs) return Term();

        // Compute the consequent at its anchor cycle (offset by any
        // `##k` on the RHS), then delay the antecedent by that
        // offset using a chain of 1-bit latches so the resulting
        // implication is expressed in the current cycle. A *range*
        // delay on a single-element consequent (`##[m:n] Q`) instead
        // becomes an OR over "Q held i cycles ago" for i spanning the
        // window, anchored at the window's latest cycle (delay +=
        // wmax) -- checking at that cycle is exactly when a violation
        // (the whole window has passed with no match) becomes certain.
        uint32_t delay =
            (b.op == BinaryAssertionOperator::NonOverlappedImplication) ? 1 : 0;
        const AssertionExpr * rhs_inner = &b.right;
        Term rhs;
        uint32_t rhs_span = 0;
        if (auto matched = match_const_delay_seq(b.right)) {
          delay += matched->first;
          rhs_inner = matched->second;
          rhs = assertion_expr_to_bool(*rhs_inner, prefix, rhs_span);
        } else if (b.right.kind == AssertionExprKind::SequenceConcat
                   && b.right.as<SequenceConcatExpr>().elements.size() == 1) {
          auto & elem = b.right.as<SequenceConcatExpr>().elements[0];
          if (!elem.delay.max) {
            // An unbounded wait is an eventuality, so there is no
            // current-cycle Boolean to return. Decline, and the
            // caller falls through to the tableau, which has an F.
            return Term();
          }
          uint32_t wmin = elem.delay.min;
          uint32_t wmax = *elem.delay.max;
          if (Term inner = assertion_expr_to_bool(*elem.sequence, prefix)) {
            delay += wmax;
            rhs = inner;
            for (uint32_t i = 1; i <= wmax - wmin; ++i) {
              rhs = solver_->make_term(
                  Or, rhs, tableau_.make_history_chain(inner, i, prefix));
            }
          }
        } else {
          rhs = assertion_expr_to_bool(*rhs_inner, prefix, rhs_span);
        }
        if (!rhs) return Term();

        // The consequent sits `delay` cycles after the attempt starts,
        // and may itself have re-anchored a further `rhs_span`; the
        // antecedent sits at `lhs_span`.  Bring both to whichever is
        // later.
        uint32_t anchor = std::max(lhs_span, delay + rhs_span);
        lhs = reanchor(lhs, lhs_span, anchor, prefix);
        rhs = reanchor(rhs, delay + rhs_span, anchor, prefix);
        Term result = solver_->make_term(Implies, lhs, rhs);
        if (lhs_delay > 0) {
          result = solver_->make_term(
              Or, tableau_.before_cycle(lhs_delay, prefix), result);
        }
        // `disable iff`: exempt this cycle's check if the disable
        // condition held anywhere in the antecedent-to-consequent
        // shift window, not just at the single cycle the check is
        // anchored at.
        if (Term dw = tableau_.disable_window(
                current_disable_cond_, anchor, prefix)) {
          result = solver_->make_term(Or, dw, result);
        }
        // Already gated above, so the caller has nothing left to do.
        span = 0;
        return result;
      }

      uint32_t lhs_span = 0;
      uint32_t rhs_span = 0;
      Term lhs = assertion_expr_to_bool(b.left, prefix, lhs_span);
      Term rhs = assertion_expr_to_bool(b.right, prefix, rhs_span);
      if (!lhs || !rhs) return Term();
      span = std::max(lhs_span, rhs_span);
      lhs = reanchor(lhs, lhs_span, span, prefix);
      rhs = reanchor(rhs, rhs_span, span, prefix);
      switch (b.op) {
        case BinaryAssertionOperator::And:
          return solver_->make_term(And, lhs, rhs);
        case BinaryAssertionOperator::Or:
          return solver_->make_term(Or, lhs, rhs);
        case BinaryAssertionOperator::Iff:
          return solver_->make_term(Equal, lhs, rhs);
        default:
          // Intersect / Throughout / Within / FollowedBy: sequence
          // operators that span multiple cycles in the general case --
          // out of scope for the current encoder. Until/SUntil/
          // UntilWith/SUntilWith aren't reducible to a current-cycle
          // Boolean either, but (unlike these) are handled by the
          // general LTL tableau in ltl_to_sat() instead.
          return Term();
      }
    }

    default:
      // FirstMatch, SequenceWithMatch, Abort, DisableIff, etc.: shapes
      // this current-cycle-Boolean fast path doesn't reduce
      // (Conditional/Case/Unary are handled by their own cases above,
      // not here). The caller throws if the ltl_to_sat() fallback
      // can't reduce it either.
      return Term();
  }
}

// ============================================================================
// Assertion-statement dispatch (formerly two StatementKind cases in
// statement.cpp)
// ============================================================================

void AssertionWalker::process_concurrent_assertion(
    const slang::ast::ConcurrentAssertionStatement & ca,
    const slang::ast::Statement & stmt,
    const string & prefix,
    const slang::ast::Expression * default_disable_expr)
{
  using namespace slang::ast;

  // Handle 'assert', 'assume', 'restrict', and 'cover'.  'assume'/
  // 'restrict' share the exact same property-shape handling as
  // 'assert' below; they differ only in what happens to the
  // resulting boolean once it's built (see the two branches further
  // down). 'cover' is handled via reachability duality (see the
  // "SVA design decisions" note at the top of this file): `cover
  // property(P)` is checked exactly like `assert property(!P)`, so
  // a "violation" of that surrogate assertion is precisely "P was
  // reached" -- it shares the same safety fast path as assert/
  // assume, just with the boolean negated before the shared
  // disable-window/push logic runs.
  // `expect (property_expr);` is a procedural blocking-wait
  // statement (pause until the property holds), not a checked
  // invariant -- a simulation-control construct with no
  // synthesizable hardware meaning, the same category as `wait`
  // elsewhere in this encoder. Handle it before the dispatch so it
  // doesn't fall through the rest of this function silently.
  if (ca.assertionKind == AssertionKind::Expect) {
    logger.log(1,
               "SystemVerilogEncoder: ignoring 'expect' property "
               "(simulation-only construct)");
    return;
  }
  bool is_assumption = ca.assertionKind == AssertionKind::Assume
                       || ca.assertionKind == AssertionKind::Restrict;
  // `cover sequence(S)` shares the exact same reachability-duality
  // treatment as `cover property(P)` below -- both just check
  // "was propertySpec ever true", regardless of whether the
  // source wrote `property` or `sequence`.
  bool is_cover = ca.assertionKind == AssertionKind::CoverProperty
                  || ca.assertionKind == AssertionKind::CoverSequence;
  if (ca.assertionKind != AssertionKind::Assert && !is_assumption
      && !is_cover) {
    return;
  }

  // Record this property's label for check_clock()'s exception message
  // -- design_clock_sym_/design_clock_edge_ themselves are NOT reset
  // here: they track the one clock the whole design is allowed to
  // have, not just this property's.
  current_assertion_label_ = assertion_label(stmt);

  // Strip the clocking wrapper (the clock event is already
  // baked into our per-cycle abstraction) and any explicit
  // `disable iff` wrapper, recording its condition.  If the
  // statement has no explicit `disable iff`, fall back to the
  // caller-resolved enclosing module's `default disable iff`, if any.
  const AssertionExpr * a = &ca.propertySpec;
  const Expression * disable_expr = nullptr;
  while (true) {
    if (a->kind == AssertionExprKind::Clocking) {
      auto & clk_expr = a->as<ClockingAssertionExpr>();
      check_clock(clk_expr.clocking);
      a = &clk_expr.expr;
    } else if (a->kind == AssertionExprKind::DisableIff) {
      auto & di = a->as<DisableIffAssertionExpr>();
      disable_expr = &di.condition;
      a = &di.expr;
    } else {
      break;
    }
  }
  if (!disable_expr) disable_expr = default_disable_expr;

  // Restored on every exit, including the throws below.  Those abort
  // the whole encode today, but one owner is one less thing to keep in
  // sync.
  struct DisableCondRestore
  {
    Term & slot;
    Term saved;
    ~DisableCondRestore() { slot = saved; }
  } disable_restore{ current_disable_cond_, current_disable_cond_ };

  if (disable_expr) {
    current_disable_cond_ = expr_encoder_.expr_to_bool(*disable_expr, prefix);
  } else {
    current_disable_cond_ = Term();
  }

  // Both encodings are built around the same thing: `violated`, the
  // condition under which the property fails at the cycle the check is
  // anchored at.  Everything decided *after* that -- cover duality,
  // the `disable iff` exemption, the per-cycle closure, which vector
  // the result lands in -- happens exactly once, below.  That is
  // deliberate: when the two encodings each applied these for
  // themselves, they twice ended up disagreeing (the LTL side dropped
  // `disable iff`, and checked only cycle 0).
  TermVec justice;
  Term violated;
  // Whether `violated` is a current-state-only predicate, which is
  // exactly assertion_expr_to_bool()'s contract and is what makes a
  // plain safety property sound.  A tableau term is not: its promise
  // inputs are pinned by a constraint one cycle later, so at the last
  // step of a bounded unrolling they are still free.
  bool per_cycle;

  // Prefer the pure-safety encoding when the property reduces to a
  // single current-cycle Boolean (plain `assert P`, `always P`,
  // bounded `|->` / `|=>` / `##k` implications).
  // assertion_expr_to_bool returns null as soon as a genuine liveness
  // operator (eventually / unbounded until) appears.
  uint32_t span = 0;
  if (Term holds = assertion_expr_to_bool(*a, prefix, span)) {
    violated = solver_->make_term(Not, holds);
    per_cycle = true;
  } else {
    span = 0;
    if (is_assumption) {
      // Temporal (non-safety) assume/restrict properties would need
      // their own fairness-constraint machinery (assuming a GF
      // condition rather than proving one), which nothing else in
      // the encoder builds yet. Dropping an assumption silently is
      // worse than dropping an assertion: the model would be left
      // *less* constrained than the source describes, so any
      // counterexample BMC/IC3 finds afterward could be spurious
      // (ruled out by the assumption this never applied) -- throw
      // rather than risk reporting an unsound "bug".
      throw PonoException(
          "SystemVerilogEncoder: temporal (non-safety) '"
          + std::string(ca.assertionKind == AssertionKind::Restrict ? "restrict"
                                                                    : "assume")
          + " property' is not supported: "
          + make_name(prefix, current_assertion_label_));
    }
    if (is_cover) {
      // The duality below would express this fine -- a cover's
      // `violated` is just ltl_to_sat(*a, /*neg=*/false, ...), and its
      // justice conditions are then the discharges of P itself, which
      // is what "P was reached" needs.  What is missing is a way to
      // tell the caller that this entry's verdict is inverted;
      // Result::ltl_justice carries no such flag.  Throw until it
      // does, rather than report a cover backwards.
      throw PonoException(
          "SystemVerilogEncoder: temporal/sequence-shaped 'cover "
          "property' is not supported");
    }
    // Build the general LTL tableau for the *negated* property and
    // collect its eventuality-discharge justice conditions.
    violated = ltl_to_sat(*a, /*neg=*/true, justice, prefix);
    if (!violated) {
      // ltl_to_sat() throws for every AssertionExprKind/operator it
      // doesn't model; a null result here can only come from a
      // bounded-sequence shape offsets_ending_now()/match_exists()
      // doesn't model (a separate, already-documented gap in that
      // primitive) -- still throw rather than silently drop.
      throw PonoException(
          "SystemVerilogEncoder: property '" + current_assertion_label_
          + "' uses an assertion shape this encoder cannot translate");
    }
    per_cycle = false;
  }

  // `cover property (P)` is checked as `assert property (!P)`, so in
  // violation space it is simply the opposite polarity.  Expressed
  // this way it no longer has to be sequenced before the exemption
  // below, the way it did when this was written against "holds".
  if (is_cover) {
    violated = solver_->make_term(Not, violated);
  }

  // If the check re-anchored itself `span` cycles forward, the first
  // `span` cycles describe attempts that could not have started, so
  // there is no verdict to give there.
  if (span > 0) {
    violated = solver_->make_term(
        And,
        solver_->make_term(Not, tableau_.before_cycle(span, prefix)),
        violated);
  }

  // `disable iff`: an attempt aborts if the condition holds anywhere
  // while it is being evaluated, so a failure there does not count.
  // How far "anywhere" reaches is what the two encodings differ on.
  if (current_disable_cond_) {
    if (per_cycle) {
      // The check was re-anchored to the last cycle it names, so the
      // attempt spans the `span` cycles ending here: exempt it if the
      // condition held at any of them.
      Term dw = tableau_.disable_window(current_disable_cond_, span, prefix);
      violated = solver_->make_term(And, solver_->make_term(Not, dw), violated);
    } else {
      // A liveness attempt is only ever violated by never completing,
      // so its evaluation runs to infinity and any later condition
      // still aborts it. Exempting the anchor cycle alone made the
      // assertion stronger than written; the attempt survives exactly
      // when the condition never holds from here on.
      violated = solver_->make_term(
          And,
          violated,
          tableau_.make_G(solver_->make_term(Not, current_disable_cond_),
                          prefix));
    }
  }

  if (is_assumption) {
    // Hold at every reachable step (init and, via the transition
    // relation, every subsequent state) -- the same "always true"
    // primitive already used for plain state/input invariants
    // elsewhere in the encoder, just applied to an assumption instead
    // of a proof obligation.
    fts_.add_constraint(solver_->make_term(Not, violated),
                        /*to_init_and_next=*/true);
    logger.log(1,
               "SystemVerilogEncoder: extracted assumption constraint "
               "from {}",
               make_name(prefix, current_assertion_label_));
    return;
  }

  if (per_cycle) {
    // A Pono safety property already means "in every reachable
    // state", so the property expression's implicit per-cycle closure
    // comes for free here.
    propvec_.push_back(solver_->make_term(Not, violated));
    logger.log(1,
               "SystemVerilogEncoder: extracted safety assertion "
               "property {} (index {})",
               make_name(prefix, current_assertion_label_),
               propvec_.size() - 1);
    return;
  }

  // Per-property activation latch: a free Boolean constant.  The
  // justice set forces it true (so this property's obligation is
  // enabled), while every *other* property's latch may stay false,
  // leaving their obligations vacuous.  This keeps independent LTL
  // properties from interfering in one system.
  Term act = fts_.make_statevar(
      make_name(prefix, "__ltl_act_" + std::to_string(tableau_.next_id())),
      solver_->make_sort(BOOL));
  fts_.assign_next(act, act);

  // The closure the safety branch got for free has to be built here:
  // the LRM evaluates a property expression at every clock tick, so a
  // counterexample is a violation *somewhere*, not only in the first
  // cycle.  Anchored at cycle 0 via the shared init flag, and added to
  // the transition relation (it references the tableau's promise
  // inputs) rather than to the initial-state predicate.
  Term obligation = solver_->make_term(
      Implies,
      solver_->make_term(And, tableau_.init_flag(prefix), act),
      tableau_.make_F(violated, justice, prefix));
  fts_.add_constraint(obligation, /*to_init_and_next=*/false);

  justice.push_back(act);
  ltl_justice_.push_back(justice);
  logger.log(1,
             "SystemVerilogEncoder: extracted LTL liveness property "
             "{} (index {}, {} justice condition(s))",
             make_name(prefix, current_assertion_label_),
             ltl_justice_.size() - 1,
             justice.size());
}

void AssertionWalker::process_immediate_assertion(
    const slang::ast::ImmediateAssertionStatement & ia,
    const smt::Term & condition,
    const string & prefix)
{
  using namespace slang::ast;

  // A *procedural* immediate assertion (`assert (expr);`),
  // distinct from the concurrent `assert property (...)` form
  // above. Reuses the same Assert-vs-Assume/Restrict split: an
  // assert becomes a safety property, an assume/restrict becomes a
  // standing constraint. Guarded by the accumulated path `condition`
  // (e.g. an enclosing `if`) rather than treated as always-active,
  // since it's only actually reached when program flow gets there --
  // "if reached, expr must hold" for assert, "if reached, assume
  // expr" for assume/restrict. Pass/fail action blocks (`assert (x)
  // else $error(...);`) are simulation-only display statements with
  // no synthesis meaning and are intentionally not processed, same as
  // $display/$error elsewhere in this encoder. `cover` uses the same
  // reachability-duality contract as process_concurrent_assertion()
  // above (see the "SVA design decisions" note at the top of this
  // file): `cover (expr);` is checked exactly like `assert (!expr);`.
  bool is_assumption = ia.assertionKind == AssertionKind::Assume
                       || ia.assertionKind == AssertionKind::Restrict;
  bool is_cover = ia.assertionKind == AssertionKind::CoverProperty;
  if (ia.assertionKind != AssertionKind::Assert && !is_assumption
      && !is_cover) {
    return;
  }
  Term bool_cond = expr_encoder_.expr_to_bool(ia.cond, prefix);
  if (is_cover) {
    bool_cond = solver_->make_term(Not, bool_cond);
  }
  Term prop = (condition == solver_->make_term(true))
                  ? bool_cond
                  : solver_->make_term(Implies, condition, bool_cond);
  if (is_assumption) {
    fts_.add_constraint(prop, /*to_init_and_next=*/true);
    logger.log(1, "SystemVerilogEncoder: extracted assumption constraint");
  } else {
    propvec_.push_back(prop);
    logger.log(1,
               "SystemVerilogEncoder: extracted safety assertion "
               "property from immediate assertion (index {})",
               propvec_.size() - 1);
  }
}

}  // namespace pono
