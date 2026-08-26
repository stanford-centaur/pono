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
 *     one `@(posedge clkN)`): this encoder has no clock-domain-crossing
 *     model at all -- every design in this frontend's test suite already
 *     implicitly assumes one global clock advances the whole design by one
 *     cycle per sample. Rather than reject a multiclock property outright,
 *     every named clock is treated identically (each clock edge mentioned
 *     anywhere in a sequence advances the same global pono-cycle). This is
 *     the minimal behavior consistent with the encoder's existing single-
 *     clock assumption, not a real multi-domain/clock-ratio model.
 */
#include "frontends/systemverilog/assertion_walker.h"

#include <cstdint>
#include <string>

#include "frontends/systemverilog/expr_encoder.h"
#include "frontends/systemverilog/tableau.h"
#include "slang/ast/Expression.h"
#include "slang/ast/expressions/AssertionExpr.h"
#include "slang/ast/expressions/MiscExpressions.h"
#include "slang/ast/statements/MiscStatements.h"
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

// A named `sequence`/`property` declaration referenced by name (e.g.
// `assert property (p_check);`) binds as a SimpleAssertionExpr wrapping
// an AssertionInstanceExpression -- not a plain boolean Expression --
// so routing it through expr_to_term() throws "unsupported expression
// kind". Slang has already expanded the referenced item's body (with
// any arguments substituted) into `body`; for the no-argument,
// non-recursive case that's exactly the AssertionExpr this encoder
// should recurse into instead. Returns nullptr if `e` isn't such a
// reference (the caller should fall back to its normal expr_to_term()
// path). Argumented, local-variable-bearing, or recursive property
// instantiations are a materially harder problem (they need their own
// binding environment) and are out of scope; throws a clear error
// rather than silently mis-evaluating them.
const slang::ast::AssertionExpr * resolve_named_assertion_ref(
    const slang::ast::Expression & e)
{
  using namespace slang::ast;
  if (e.kind != ExpressionKind::AssertionInstance) return nullptr;
  auto & aie = e.as<AssertionInstanceExpression>();
  if (aie.isRecursiveProperty || !aie.arguments.empty()
      || !aie.localVars.empty()) {
    throw PonoException(
        "SystemVerilogEncoder: named sequence/property references with "
        "arguments, local variables, or recursion are not supported");
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
    const slang::ast::AssertionExpr & seq, const string & prefix)
{
  using namespace slang::ast;

  // A single Boolean expression, optionally with a consecutive
  // repetition (`expr[*n:m]`, `expr[+]`, `expr[*]`). Shared by both
  // SimpleAssertionExpr and SequenceWithMatchExpr, which each carry
  // their own std::optional<SequenceRepetition>.
  auto boolean_with_repetition =
      [&](const slang::ast::Expression & expr,
          const std::optional<SequenceRepetition> & repetition) -> TermVec {
    Term b = expr_encoder_.expr_to_term(expr, prefix);
    b = solver_->make_term(Distinct, b, solver_->make_term(0, b->get_sort()));
    if (!repetition) return { b };
    if (repetition->kind != SequenceRepetition::Consecutive) {
      throw PonoException(
          "SystemVerilogEncoder: nonconsecutive/goto sequence repetition "
          "([->]/[=]) is not supported");
    }
    if (!repetition->range.max) {
      throw PonoException(
          "SystemVerilogEncoder: unbounded consecutive sequence "
          "repetition ([*]/[+]/[*n:$]) is not supported");
    }
    uint32_t lo = repetition->range.min;
    uint32_t hi = *repetition->range.max;
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
            And, running, tableau_.delay_bool(b, count - 1, prefix));
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
      result =
          solver_->make_term(Or, result, tableau_.delay_bool(base, j, prefix));
    }
    return result;
  };

  // AND of `base` and its delayed copies over the last `k` cycles --
  // "did `base` hold at *every* point in the last k+1 cycles". Used by
  // Throughout.
  auto window_and = [&](const Term & base, uint32_t k) -> Term {
    Term result = base;
    for (uint32_t j = 1; j <= k; ++j) {
      result =
          solver_->make_term(And, result, tableau_.delay_bool(base, j, prefix));
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
      return offsets_ending_now(swm.expr, prefix);
    }

    case AssertionExprKind::FirstMatch:
      // first_match(seq) only restricts *which* match is reported when
      // a sequence can match in more than one way -- it never changes
      // whether a match exists at all, which is all this encoder's
      // callers (an implication antecedent, an intersect/within/
      // throughout operand) ever ask offsets_ending_now() for.
      return offsets_ending_now(seq.as<FirstMatchAssertionExpr>().seq, prefix);

    case AssertionExprKind::Clocking:
      // Per this file's multiclock design decision: every named clock
      // is treated as the same global pono-cycle, so a nested clocking
      // change inside a sequence element is simply unwrapped.
      return offsets_ending_now(seq.as<ClockingAssertionExpr>().expr, prefix);

    case AssertionExprKind::SequenceConcat: {
      auto & sc = seq.as<SequenceConcatExpr>();
      TermVec acc;
      for (size_t i = 0; i < sc.elements.size(); ++i) {
        auto & elem = sc.elements[i];
        if (!elem.delay.max) {
          throw PonoException(
              "SystemVerilogEncoder: unbounded sequence delay (##[m:$]) is "
              "not supported");
        }
        uint32_t dmin = elem.delay.min;
        uint32_t dmax = *elem.delay.max;
        TermVec elem_offsets = offsets_ending_now(*elem.sequence, prefix);
        if (elem_offsets.empty()) return {};

        if (i == 0) {
          // The delay before the very first element just relabels how
          // far back "the sequence's start" is, with no extra
          // condition to AND in.
          size_t new_size = elem_offsets.size() - 1 + dmax + 1;
          if (new_size > MAX_SEQ_WINDOW) {
            throw PonoException("SystemVerilogEncoder: sequence window exceeds "
                                + std::to_string(MAX_SEQ_WINDOW) + " cycles");
          }
          acc.assign(new_size, Term());
          for (size_t l = 0; l < elem_offsets.size(); ++l) {
            if (!elem_offsets[l]) continue;
            for (uint32_t d = dmin; d <= dmax; ++d) {
              size_t idx = l + d;
              acc[idx] = acc[idx]
                             ? solver_->make_term(Or, acc[idx], elem_offsets[l])
                             : elem_offsets[l];
            }
          }
          continue;
        }

        size_t new_size = acc.size() - 1 + dmax + (elem_offsets.size() - 1) + 1;
        if (new_size > MAX_SEQ_WINDOW) {
          throw PonoException("SystemVerilogEncoder: sequence window exceeds "
                              + std::to_string(MAX_SEQ_WINDOW) + " cycles");
        }
        TermVec new_acc(new_size, Term());
        for (size_t lp = 0; lp < acc.size(); ++lp) {
          if (!acc[lp]) continue;
          for (uint32_t d = dmin; d <= dmax; ++d) {
            for (size_t le = 0; le < elem_offsets.size(); ++le) {
              if (!elem_offsets[le]) continue;
              size_t idx = lp + d + le;
              // Bring the (already-anchored-at-"now") prefix condition
              // back by (d + le) cycles so it aligns with this
              // element's own occurrence, then AND with this
              // element's own (unshifted) completion condition.
              Term shifted_prefix =
                  (d + le == 0) ? acc[lp]
                                : tableau_.delay_bool(acc[lp], d + le, prefix);
              Term combined =
                  solver_->make_term(And, shifted_prefix, elem_offsets[le]);
              new_acc[idx] =
                  new_acc[idx] ? solver_->make_term(Or, new_acc[idx], combined)
                               : combined;
            }
          }
        }
        acc = std::move(new_acc);
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
        throw PonoException(
            "SystemVerilogEncoder: weak()/strong() of a sequence with its "
            "own leading repetition is not supported");
      }
      Term t = expr_encoder_.expr_to_term(simple.expr, prefix);
      return solver_->make_term(
          Distinct, t, solver_->make_term(0, t->get_sort()));
    }
    case AssertionExprKind::FirstMatch:
      return leading_condition(seq.as<FirstMatchAssertionExpr>().seq, prefix);
    case AssertionExprKind::Clocking:
      return leading_condition(seq.as<ClockingAssertionExpr>().expr, prefix);
    case AssertionExprKind::SequenceConcat:
      return leading_condition(
          *seq.as<SequenceConcatExpr>().elements[0].sequence, prefix);
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
      tableau_.delay_bool(leading_condition(seq, prefix), s, prefix);
  Term completed_in_window = me;
  for (uint32_t j = 1; j <= s; ++j) {
    completed_in_window = solver_->make_term(
        Or, completed_in_window, tableau_.delay_bool(me, j, prefix));
  }
  // Violated iff an attempt began exactly S cycles ago and no
  // completion happened anywhere from then through now; weak(seq) is
  // the negation -- no obligation to ever attempt, but an attempt that
  // did begin must not be a definite, provable failure.
  Term violated = solver_->make_term(
      And, started_s_ago, solver_->make_term(Not, completed_in_window));
  return solver_->make_term(Not, violated);
}

smt::Term AssertionWalker::ltl_to_sat(const slang::ast::AssertionExpr & ae,
                                      bool neg,
                                      smt::TermVec & justice,
                                      const string & prefix)
{
  using namespace slang::ast;

  switch (ae.kind) {
    case AssertionExprKind::Clocking:
      return ltl_to_sat(
          ae.as<ClockingAssertionExpr>().expr, neg, justice, prefix);

    case AssertionExprKind::StrongWeak: {
      auto & sw = ae.as<StrongWeakAssertionExpr>();
      // The strong/weak qualifier only meaningfully differs for a
      // genuine bounded sequence -- must it eventually complete
      // (strong) or not (weak)? Any other shape (already a plain
      // Boolean/temporal expression) is unaffected by the qualifier
      // under this encoder's infinite-lasso semantics; just unwrap.
      if (sw.strength == StrongWeakAssertionExpr::Strong) {
        Term me = match_exists(sw.expr, prefix);
        if (me) {
          // strong(seq): a genuine liveness obligation -- the sequence
          // must eventually complete a match.
          return neg ? tableau_.make_G(solver_->make_term(Not, me), prefix)
                     : tableau_.make_F(me, justice, prefix);
        }
      }
      return ltl_to_sat(sw.expr, neg, justice, prefix);
    }

    case AssertionExprKind::Simple: {
      auto & simple = ae.as<SimpleAssertionExpr>();
      if (auto * named = resolve_named_assertion_ref(simple.expr)) {
        return ltl_to_sat(*named, neg, justice, prefix);
      }
      if (simple.repetition) {
        // See the matching check in assertion_expr_to_bool(): route
        // through the general bounded sequence matcher (which throws
        // for an unbounded repeat count) instead of silently ignoring
        // the repetition.
        Term me = match_exists(ae, prefix);
        if (!me) return Term();
        return neg ? solver_->make_term(Not, me) : me;
      }
      Term t = expr_encoder_.expr_to_term(simple.expr, prefix);
      if (!t) return Term();
      Term zero = solver_->make_term(0, t->get_sort());
      Term b = solver_->make_term(Distinct, t, zero);
      return neg ? solver_->make_term(Not, b) : b;
    }

    case AssertionExprKind::SequenceConcat: {
      // A bare `##k Q` property has the same truth value as Q under
      // our infinite-time semantics (modulo a front shift).  Unwrap.
      if (auto m = match_const_delay_seq(ae)) {
        return ltl_to_sat(*m->second, neg, justice, prefix);
      }
      return Term();
    }

    case AssertionExprKind::Unary: {
      auto & u = ae.as<UnaryAssertionExpr>();
      switch (u.op) {
        case UnaryAssertionOperator::Not:
          return ltl_to_sat(u.expr, !neg, justice, prefix);

        case UnaryAssertionOperator::Always:
        case UnaryAssertionOperator::SAlways: {
          // G phi  (positive)  /  !G phi == F !phi  (negated)
          Term phi = ltl_to_sat(u.expr, neg, justice, prefix);
          if (!phi) return Term();
          return neg ? tableau_.make_F(phi, justice, prefix)
                     : tableau_.make_G(phi, prefix);
        }

        case UnaryAssertionOperator::Eventually:
        case UnaryAssertionOperator::SEventually: {
          // F phi  (positive)  /  !F phi == G !phi  (negated)
          Term phi = ltl_to_sat(u.expr, neg, justice, prefix);
          if (!phi) return Term();
          return neg ? tableau_.make_G(phi, prefix)
                     : tableau_.make_F(phi, justice, prefix);
        }

        case UnaryAssertionOperator::NextTime:
        case UnaryAssertionOperator::SNextTime: {
          // !X phi == X !phi, so the negation rides along inside phi.
          Term phi = ltl_to_sat(u.expr, neg, justice, prefix);
          if (!phi) return Term();
          return tableau_.make_X(phi, prefix);
        }

        default: return Term();
      }
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
          if (auto m = match_const_delay_seq(b.right)) {
            delay += m->first;
            rhs = m->second;
          }
          Term l = ltl_to_sat(b.left, !neg, justice, prefix);
          Term r = ltl_to_sat(*rhs, neg, justice, prefix);
          if (!l || !r) return Term();
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

        default:
          // Intersect / Throughout / Within / FollowedBy: multi-cycle
          // sequence operators the tableau does not model.
          return Term();
      }
    }

    default: return Term();
  }
}

smt::Term AssertionWalker::assertion_expr_to_bool(
    const slang::ast::AssertionExpr & ae, const string & prefix)
{
  using namespace slang::ast;

  switch (ae.kind) {
    case AssertionExprKind::Clocking: {
      // The clocking event has already been baked into our cycle
      // abstraction; just recurse into the underlying expression.
      return assertion_expr_to_bool(ae.as<ClockingAssertionExpr>().expr,
                                    prefix);
    }

    case AssertionExprKind::Simple: {
      auto & simple = ae.as<SimpleAssertionExpr>();
      if (auto * named = resolve_named_assertion_ref(simple.expr)) {
        return assertion_expr_to_bool(*named, prefix);
      }
      if (simple.repetition) {
        // `expr[*n:m]`/`expr[+]`/`expr[*]`: route through the general
        // bounded sequence matcher (which throws for an unbounded
        // repeat count) instead of silently ignoring the repetition
        // and returning a bare `bool(expr)`.
        return match_exists(ae, prefix);
      }
      Term t = expr_encoder_.expr_to_term(simple.expr, prefix);
      // Normalize to Bool: t != 0.
      Sort sort = t->get_sort();
      Term zero = solver_->make_term(0, sort);
      return solver_->make_term(Distinct, t, zero);
    }

    case AssertionExprKind::SequenceConcat: {
      // A standalone `##k Q` as a property: in our infinite-time
      // safety encoding the constant front-shift doesn't change the
      // truth value (it just postpones when the first violation can
      // be reported), so unwrap the inner sequence.
      if (auto matched = match_const_delay_seq(ae)) {
        return assertion_expr_to_bool(*matched->second, prefix);
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
      switch (u.op) {
        case UnaryAssertionOperator::Not: {
          Term inner = assertion_expr_to_bool(u.expr, prefix);
          if (!inner) return Term();
          return solver_->make_term(Not, inner);
        }
        case UnaryAssertionOperator::Always:
        case UnaryAssertionOperator::SAlways: {
          // Pure safety: `always P` is true at the current cycle
          // exactly when P is true at the current cycle (the
          // "always at every cycle" closure is implicit in the
          // per-cycle property check).
          return assertion_expr_to_bool(u.expr, prefix);
        }
        default:
          // Eventually / SEventually / NextTime / SNextTime can't
          // be folded into a current-cycle Boolean: liveness must
          // come through the top-level dispatch and NextTime needs
          // a forward-shift the encoder doesn't model yet.
          return Term();
      }
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
        Term lhs = assertion_expr_to_bool(*lhs_inner, prefix);
        if (!lhs) lhs = match_exists(*lhs_inner, prefix);
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
        if (auto matched = match_const_delay_seq(b.right)) {
          delay += matched->first;
          rhs_inner = matched->second;
          rhs = assertion_expr_to_bool(*rhs_inner, prefix);
        } else if (b.right.kind == AssertionExprKind::SequenceConcat
                   && b.right.as<SequenceConcatExpr>().elements.size() == 1) {
          auto & elem = b.right.as<SequenceConcatExpr>().elements[0];
          if (!elem.delay.max) {
            throw PonoException(
                "SystemVerilogEncoder: unbounded sequence delay (##[m:$]) "
                "is not supported");
          }
          uint32_t wmin = elem.delay.min;
          uint32_t wmax = *elem.delay.max;
          if (Term inner = assertion_expr_to_bool(*elem.sequence, prefix)) {
            delay += wmax;
            rhs = inner;
            for (uint32_t i = 1; i <= wmax - wmin; ++i) {
              rhs = solver_->make_term(
                  Or, rhs, tableau_.delay_bool(inner, i, prefix));
            }
          }
        } else {
          rhs = assertion_expr_to_bool(*rhs_inner, prefix);
        }
        if (!rhs) return Term();

        if (delay > 0) {
          // Materialize the antecedent as a 1-bit BV so the
          // history chain has a value to latch.
          Sort bv1 = solver_->make_sort(BV, 1);
          Term one_bv1 = solver_->make_term(1, bv1);
          Term zero_bv1 = solver_->make_term(0, bv1);
          Term lhs_bv = solver_->make_term(Ite, lhs, one_bv1, zero_bv1);
          Term delayed_bv = tableau_.make_history_chain(lhs_bv, delay, prefix);
          lhs = solver_->make_term(Equal, delayed_bv, one_bv1);
        }
        Term result = solver_->make_term(Implies, lhs, rhs);
        if (lhs_delay > 0) {
          result = solver_->make_term(
              Or, tableau_.before_cycle(lhs_delay, prefix), result);
        }
        // `disable iff`: exempt this cycle's check if the disable
        // condition held anywhere in the antecedent-to-consequent
        // shift window, not just at the single cycle the check is
        // anchored at.
        if (Term dw =
                tableau_.disable_window(current_disable_cond_, delay, prefix)) {
          result = solver_->make_term(Or, dw, result);
        }
        return result;
      }

      Term lhs = assertion_expr_to_bool(b.left, prefix);
      Term rhs = assertion_expr_to_bool(b.right, prefix);
      if (!lhs || !rhs) return Term();
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
      // FirstMatch, SequenceWithMatch, Conditional, Case, Abort,
      // DisableIff, etc.: shapes this current-cycle-Boolean fast path
      // doesn't reduce (Unary is handled by its own case above, not
      // here).  Caller logs the skipped kind if the ltl_to_sat()
      // fallback can't reduce it either.
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

  // Strip the clocking wrapper (the clock event is already
  // baked into our per-cycle abstraction) and any explicit
  // `disable iff` wrapper, recording its condition.  If the
  // statement has no explicit `disable iff`, fall back to the
  // caller-resolved enclosing module's `default disable iff`, if any.
  const AssertionExpr * a = &ca.propertySpec;
  const Expression * disable_expr = nullptr;
  while (true) {
    if (a->kind == AssertionExprKind::Clocking) {
      a = &a->as<ClockingAssertionExpr>().expr;
    } else if (a->kind == AssertionExprKind::DisableIff) {
      auto & di = a->as<DisableIffAssertionExpr>();
      disable_expr = &di.condition;
      a = &di.expr;
    } else {
      break;
    }
  }
  if (!disable_expr) disable_expr = default_disable_expr;

  Term saved_disable_cond = current_disable_cond_;
  if (disable_expr) {
    Term dc = expr_encoder_.expr_to_term(*disable_expr, prefix);
    current_disable_cond_ =
        solver_->make_term(Distinct, dc, solver_->make_term(0, dc->get_sort()));
  } else {
    current_disable_cond_ = Term();
  }

  // Prefer the pure-safety encoding when the property reduces to
  // a single current-cycle Boolean (plain `assert P`, `always
  // P`, bounded `|->` / `|=>` / `##k` implications).
  // assertion_expr_to_bool returns null as soon as a genuine
  // liveness operator (eventually / unbounded until) appears.
  if (Term prop = assertion_expr_to_bool(*a, prefix)) {
    if (is_cover) {
      // Reachability duality: negate before the disable-window
      // exemption below runs, so a `disable iff C` on a cover
      // property correctly means "don't count P as covered while
      // C holds" (the composed surrogate is `!P || C`, whose own
      // violation is exactly "P held and C did not").
      prop = solver_->make_term(Not, prop);
    }
    // Catch-all `disable iff` exemption for property shapes
    // (plain `assert P`, `always P`, ...) that don't already
    // gate themselves more precisely inside assertion_expr_to_bool.
    // For an assumption, this is exactly the right shape too:
    // "assume P disable iff C" means P is only assumed while C is
    // false, i.e. the ever-true constraint is (C || P).
    if (Term dw = tableau_.disable_window(current_disable_cond_, 0, prefix)) {
      prop = solver_->make_term(Or, dw, prop);
    }
    if (is_assumption) {
      // Hold at every reachable step (init and, via the transition
      // relation, every subsequent state) -- the same "always true"
      // primitive already used for plain state/input invariants
      // elsewhere in the encoder, just applied to an assumption
      // instead of a proof obligation.
      fts_.add_constraint(prop, /*to_init_and_next=*/true);
      logger.log(1,
                 "SystemVerilogEncoder: extracted assumption "
                 "constraint from {}",
                 make_name(prefix, assertion_label(stmt)));
    } else {
      propvec_.push_back(prop);
      logger.log(1,
                 "SystemVerilogEncoder: extracted safety assertion "
                 "property {} (index {})",
                 make_name(prefix, assertion_label(stmt)),
                 propvec_.size() - 1);
    }
    current_disable_cond_ = saved_disable_cond;
    return;
  }

  if (is_assumption) {
    // Temporal (non-safety) assume/restrict properties would need
    // their own fairness-constraint machinery (assuming a GF
    // condition rather than proving one), which nothing else in
    // the encoder builds yet -- skip cleanly rather than attempt a
    // partial, likely-unsound translation.
    logger.log(1,
               "SystemVerilogEncoder: skipping unsupported temporal "
               "assume/restrict property {}",
               make_name(prefix, assertion_label(stmt)));
    current_disable_cond_ = saved_disable_cond;
    return;
  }

  if (is_cover) {
    // A temporal/sequence-shaped cover goal (e.g. `cover property
    // (a ##1 b)`) would need the reachability duality above
    // extended through the same LTL tableau `ltl_to_sat()` builds
    // for `assert`/`assume` -- negating a liveness obligation
    // doesn't correspond to "was this ever reached" the way it
    // does for a plain current-cycle Boolean, so that extension
    // is out of scope here. Throw rather than silently drop.
    throw PonoException(
        "SystemVerilogEncoder: temporal/sequence-shaped 'cover "
        "property' is not supported");
  }

  // Otherwise build the general LTL tableau for the *negated*
  // property and collect its eventuality-discharge justice
  // conditions.  A fair lasso of the resulting system (every
  // justice condition true infinitely often) on which the
  // negated property holds at cycle 0 is exactly a
  // counterexample to the original assertion.
  TermVec justice;
  Term satpsi = ltl_to_sat(*a, /*neg=*/true, justice, prefix);
  if (!satpsi) {
    logger.log(1,
               "SystemVerilogEncoder: skipping unsupported temporal "
               "assertion kind {}",
               static_cast<int>(a->kind));
    current_disable_cond_ = saved_disable_cond;
    return;
  }

  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);

  // Per-property activation latch: a free 1-bit constant.  The
  // justice set forces it to 1 (so this property's time-0
  // obligation is enabled), while every *other* property's latch
  // may stay 0, leaving their obligations vacuous.  This keeps
  // independent LTL properties from interfering in one system.
  Term act = fts_.make_statevar(
      make_name(prefix, "__ltl_act_" + std::to_string(tableau_.next_id())),
      bv1);
  fts_.assign_next(act, act);
  Term act_bool = solver_->make_term(Equal, act, one_bv1);

  // Time-0 obligation: when active, the negated property must
  // hold at the first cycle.  Gated by the shared init flag so
  // it constrains only cycle 0, and added to the transition
  // relation (it references the tableau's promise inputs) rather
  // than to the initial-state predicate.
  Term obligation = solver_->make_term(
      Implies,
      solver_->make_term(And, tableau_.init_flag(prefix), act_bool),
      satpsi);
  fts_.add_constraint(obligation, /*to_init_and_next=*/false);

  justice.push_back(act_bool);
  ltl_justice_.push_back(justice);
  logger.log(1,
             "SystemVerilogEncoder: extracted LTL liveness property "
             "{} (index {}, {} justice condition(s))",
             make_name(prefix, assertion_label(stmt)),
             ltl_justice_.size() - 1,
             justice.size());
  current_disable_cond_ = saved_disable_cond;
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
  Term cond_term = expr_encoder_.expr_to_term(ia.cond, prefix);
  Term bool_cond = solver_->make_term(
      Distinct, cond_term, solver_->make_term(0, cond_term->get_sort()));
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
