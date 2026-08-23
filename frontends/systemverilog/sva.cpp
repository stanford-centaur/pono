/*! \file sva.cpp
 *  \brief SystemVerilogEncoder's SVA/LTL encoding: the history-chain and
 *         `$past`-style helpers, the bounded LTL-to-SAT tableau
 *         (X/G/F/R/U operators, sequence matching), and the top-level
 *         assertion-expression-to-bool conversion.
 *
 * SVA design decisions
 * ---------------------
 * A few SVA constructs have no single "obviously correct" encoding given
 * this encoder's model (a single global clock, and a propvec()-of-safety-
 * properties / ltl_justice()-of-liveness-obligations interface with no
 * notion of coverage or multiple clock domains). The choices made here are
 * deliberate and documented so future changes don't accidentally drift
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
#include <cstdint>
#include <string>

#include "frontends/systemverilog/encoder.h"
#include "slang/ast/Expression.h"
#include "slang/ast/Symbol.h"
#include "slang/ast/expressions/AssertionExpr.h"
#include "slang/ast/expressions/MiscExpressions.h"
#include "smt-switch/smt.h"
#include "utils/exceptions.h"
#include "utils/logger.h"

using namespace smt;
using namespace std;

namespace pono {

Term SystemVerilogEncoder::make_history_chain(const Term & value, uint32_t n)
{
  Sort sort = value->get_sort();
  Term zero = solver_->make_term(0, sort);
  Term link = value;
  for (uint32_t i = 0; i < n; ++i) {
    Term latch = fts_.make_statevar(
        make_name("__sva_past_" + std::to_string(latch_counter_++)), sort);
    fts_.constrain_init(solver_->make_term(Equal, latch, zero));
    fts_.assign_next(latch, link);
    link = latch;
  }
  return link;
}

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
}  // namespace

// ============================================================================
// LTL tableau
// ============================================================================
//
// Properties that are not pure safety are translated with a standard
// symbolic LTL tableau (temporal testers).  Every temporal operator
// introduces a one-step "promise" tester:
//
//   * a free 1-bit input  n  guessing whether the operator's body
//     holds at the *next* cycle (this is the operator's value now);
//   * a 1-bit state  s  with  s' = n  (it remembers the guess), plus
//     the per-cycle consistency constraint  s == body, which forces
//     the guess made one cycle earlier to have been correct.
//
// Greatest-fixpoint operators (G, R) need no fairness.  Least-fixpoint
// operators (F, strong-U) additionally emit a justice (GF) condition
// that discharges the eventuality, ruling out lassos where a promise
// is made forever but never fulfilled.  ltl_to_sat() pushes negation
// to the leaves on the fly, so the testers it builds always match the
// negation-normal form of the (negated) property.

smt::Term SystemVerilogEncoder::ltl_init_flag()
{
  if (ltl_init_flag_) return ltl_init_flag_;
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  Term flag = fts_.make_statevar(make_name("__ltl_init_flag"), bv1);
  fts_.constrain_init(solver_->make_term(Equal, flag, one_bv1));
  fts_.assign_next(flag, zero_bv1);
  ltl_init_flag_ = solver_->make_term(Equal, flag, one_bv1);
  return ltl_init_flag_;
}

smt::Term SystemVerilogEncoder::ltl_before_cycle(uint32_t k)
{
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  // "is cycle 0" pulse, delayed by 1..k-1 steps gives "is cycle i" for
  // each i in [1, k); their disjunction (with the undelayed pulse) is
  // true exactly during cycles [0, k).
  Term pulse = ltl_init_flag();
  Term result = pulse;
  Term link = solver_->make_term(Ite, pulse, one_bv1, zero_bv1);
  for (uint32_t i = 1; i < k; ++i) {
    Term latch = fts_.make_statevar(
        make_name("__before_cycle_" + std::to_string(latch_counter_++)), bv1);
    fts_.constrain_init(solver_->make_term(Equal, latch, zero_bv1));
    fts_.assign_next(latch, link);
    result = solver_->make_term(
        Or, result, solver_->make_term(Equal, latch, one_bv1));
    link = latch;
  }
  return result;
}

smt::Term SystemVerilogEncoder::disable_window(uint32_t window)
{
  if (!current_disable_cond_) return Term();
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  // OR the disable condition together with its 1..window-cycle-delayed
  // versions, so a `disable iff` that was true anywhere in an
  // antecedent's shift window still exempts the whole property, not
  // just the single cycle the check ends up anchored at.
  Term result = current_disable_cond_;
  Term link = solver_->make_term(Ite, current_disable_cond_, one_bv1, zero_bv1);
  for (uint32_t i = 0; i < window; ++i) {
    Term latch = fts_.make_statevar(
        make_name("__disable_hist_" + std::to_string(latch_counter_++)), bv1);
    fts_.constrain_init(solver_->make_term(Equal, latch, zero_bv1));
    fts_.assign_next(latch, link);
    result = solver_->make_term(
        Or, result, solver_->make_term(Equal, latch, one_bv1));
    link = latch;
  }
  return result;
}

smt::Term SystemVerilogEncoder::ltl_make_X(const smt::Term & phi)
{
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  uint32_t id = latch_counter_++;
  Term n = fts_.make_inputvar(make_name("__ltl_n_" + std::to_string(id)), bv1);
  Term s = fts_.make_statevar(make_name("__ltl_s_" + std::to_string(id)), bv1);
  fts_.assign_next(s, n);  // s@(t+1) = n@t  (remember the guess)
  // The guess made at t-1 about phi@t must have been correct.  The
  // cycle-0 instance only pins the otherwise-unused initial s.
  Term phi_bv = solver_->make_term(Ite, phi, one_bv1, zero_bv1);
  fts_.add_constraint(solver_->make_term(Equal, s, phi_bv),
                      /*to_init_and_next=*/false);
  return solver_->make_term(Equal, n, one_bv1);  // sat(X phi) = guess
}

smt::Term SystemVerilogEncoder::ltl_make_G(const smt::Term & phi)
{
  // G phi == phi && X(G phi)
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  uint32_t id = latch_counter_++;
  Term n = fts_.make_inputvar(make_name("__ltl_n_" + std::to_string(id)), bv1);
  Term s = fts_.make_statevar(make_name("__ltl_s_" + std::to_string(id)), bv1);
  fts_.assign_next(s, n);
  Term n_bool = solver_->make_term(Equal, n, one_bv1);
  Term body = solver_->make_term(And, phi, n_bool);  // sat(G phi)
  Term body_bv = solver_->make_term(Ite, body, one_bv1, zero_bv1);
  fts_.add_constraint(solver_->make_term(Equal, s, body_bv),
                      /*to_init_and_next=*/false);
  return body;
}

smt::Term SystemVerilogEncoder::ltl_make_F(const smt::Term & phi,
                                           smt::TermVec & justice)
{
  // F phi == phi || X(F phi)
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  uint32_t id = latch_counter_++;
  Term n = fts_.make_inputvar(make_name("__ltl_n_" + std::to_string(id)), bv1);
  Term s = fts_.make_statevar(make_name("__ltl_s_" + std::to_string(id)), bv1);
  fts_.assign_next(s, n);
  Term n_bool = solver_->make_term(Equal, n, one_bv1);
  Term body = solver_->make_term(Or, phi, n_bool);  // sat(F phi)
  Term body_bv = solver_->make_term(Ite, body, one_bv1, zero_bv1);
  fts_.add_constraint(solver_->make_term(Equal, s, body_bv),
                      /*to_init_and_next=*/false);
  // Discharge: infinitely often phi holds or we stop promising F.
  justice.push_back(
      solver_->make_term(Or, phi, solver_->make_term(Not, n_bool)));
  return body;
}

smt::Term SystemVerilogEncoder::ltl_make_R(const smt::Term & a,
                                           const smt::Term & b)
{
  // a R b == b && (a || X(a R b))   (greatest fixpoint, no fairness)
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  uint32_t id = latch_counter_++;
  Term n = fts_.make_inputvar(make_name("__ltl_n_" + std::to_string(id)), bv1);
  Term s = fts_.make_statevar(make_name("__ltl_s_" + std::to_string(id)), bv1);
  fts_.assign_next(s, n);
  Term n_bool = solver_->make_term(Equal, n, one_bv1);
  Term body = solver_->make_term(
      And, b, solver_->make_term(Or, a, n_bool));  // sat(a R b)
  Term body_bv = solver_->make_term(Ite, body, one_bv1, zero_bv1);
  fts_.add_constraint(solver_->make_term(Equal, s, body_bv),
                      /*to_init_and_next=*/false);
  return body;
}

smt::Term SystemVerilogEncoder::ltl_make_U(const smt::Term & a,
                                           const smt::Term & b,
                                           smt::TermVec & justice)
{
  // a U b (strong) == b || (a && X(a U b))
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  uint32_t id = latch_counter_++;
  Term n = fts_.make_inputvar(make_name("__ltl_n_" + std::to_string(id)), bv1);
  Term s = fts_.make_statevar(make_name("__ltl_s_" + std::to_string(id)), bv1);
  fts_.assign_next(s, n);
  Term n_bool = solver_->make_term(Equal, n, one_bv1);
  Term body = solver_->make_term(
      Or, b, solver_->make_term(And, a, n_bool));  // sat(a U b)
  Term body_bv = solver_->make_term(Ite, body, one_bv1, zero_bv1);
  fts_.add_constraint(solver_->make_term(Equal, s, body_bv),
                      /*to_init_and_next=*/false);
  // Discharge: infinitely often b holds or a U b is false.
  justice.push_back(solver_->make_term(Or, solver_->make_term(Not, body), b));
  return body;
}

Term SystemVerilogEncoder::delay_bool(const Term & cond, uint32_t n)
{
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  Term cond_bv = solver_->make_term(Ite, cond, one_bv1, zero_bv1);
  Term delayed_bv = make_history_chain(cond_bv, n);
  return solver_->make_term(Equal, delayed_bv, one_bv1);
}

namespace {
// The largest total span (in cycles) offsets_ending_now() will build
// before giving up -- a defensive cap against a pathological/absurd
// bounded sequence, mirroring the MAX_ITERS-style caps used elsewhere
// in this encoder for other compile-time-unrolled constructs.
constexpr uint32_t MAX_SEQ_WINDOW = 256;
}  // namespace

smt::TermVec SystemVerilogEncoder::offsets_ending_now(
    const slang::ast::AssertionExpr & seq)
{
  using namespace slang::ast;

  // A single Boolean expression, optionally with a consecutive
  // repetition (`expr[*n:m]`, `expr[+]`, `expr[*]`). Shared by both
  // SimpleAssertionExpr and SequenceWithMatchExpr, which each carry
  // their own std::optional<SequenceRepetition>.
  auto boolean_with_repetition =
      [&](const slang::ast::Expression & expr,
          const std::optional<SequenceRepetition> & repetition) -> TermVec {
    Term b = expr_to_term(expr);
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
        running = solver_->make_term(And, running, delay_bool(b, count - 1));
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
      result = solver_->make_term(Or, result, delay_bool(base, j));
    }
    return result;
  };

  // AND of `base` and its delayed copies over the last `k` cycles --
  // "did `base` hold at *every* point in the last k+1 cycles". Used by
  // Throughout.
  auto window_and = [&](const Term & base, uint32_t k) -> Term {
    Term result = base;
    for (uint32_t j = 1; j <= k; ++j) {
      result = solver_->make_term(And, result, delay_bool(base, j));
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
      return offsets_ending_now(swm.expr);
    }

    case AssertionExprKind::FirstMatch:
      // first_match(seq) only restricts *which* match is reported when
      // a sequence can match in more than one way -- it never changes
      // whether a match exists at all, which is all this encoder's
      // callers (an implication antecedent, an intersect/within/
      // throughout operand) ever ask offsets_ending_now() for.
      return offsets_ending_now(seq.as<FirstMatchAssertionExpr>().seq);

    case AssertionExprKind::Clocking:
      // Per this file's multiclock design decision: every named clock
      // is treated as the same global pono-cycle, so a nested clocking
      // change inside a sequence element is simply unwrapped.
      return offsets_ending_now(seq.as<ClockingAssertionExpr>().expr);

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
        TermVec elem_offsets = offsets_ending_now(*elem.sequence);
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
                  (d + le == 0) ? acc[lp] : delay_bool(acc[lp], d + le);
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
          TermVec v1 = offsets_ending_now(b.left);
          TermVec v2 = offsets_ending_now(b.right);
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
          TermVec v1 = offsets_ending_now(b.left);
          TermVec v2 = offsets_ending_now(b.right);
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
          Term expr_bool = assertion_expr_to_bool(b.left);
          if (!expr_bool) return {};
          TermVec v2 = offsets_ending_now(b.right);
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

smt::Term SystemVerilogEncoder::match_exists(
    const slang::ast::AssertionExpr & seq)
{
  TermVec offsets = offsets_ending_now(seq);
  Term result;
  for (auto & t : offsets) {
    if (!t) continue;
    result = result ? solver_->make_term(Or, result, t) : t;
  }
  return result;
}

smt::Term SystemVerilogEncoder::leading_condition(
    const slang::ast::AssertionExpr & seq)
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
      Term t = expr_to_term(simple.expr);
      return solver_->make_term(
          Distinct, t, solver_->make_term(0, t->get_sort()));
    }
    case AssertionExprKind::FirstMatch:
      return leading_condition(seq.as<FirstMatchAssertionExpr>().seq);
    case AssertionExprKind::Clocking:
      return leading_condition(seq.as<ClockingAssertionExpr>().expr);
    case AssertionExprKind::SequenceConcat:
      return leading_condition(
          *seq.as<SequenceConcatExpr>().elements[0].sequence);
    default:
      throw PonoException(
          "SystemVerilogEncoder: weak()/strong() of this sequence shape is "
          "not supported");
  }
}

smt::Term SystemVerilogEncoder::weak_seq_bool(
    const slang::ast::AssertionExpr & seq)
{
  TermVec offsets = offsets_ending_now(seq);
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
  Term started_s_ago = delay_bool(leading_condition(seq), s);
  Term completed_in_window = me;
  for (uint32_t j = 1; j <= s; ++j) {
    completed_in_window =
        solver_->make_term(Or, completed_in_window, delay_bool(me, j));
  }
  // Violated iff an attempt began exactly S cycles ago and no
  // completion happened anywhere from then through now; weak(seq) is
  // the negation -- no obligation to ever attempt, but an attempt that
  // did begin must not be a definite, provable failure.
  Term violated = solver_->make_term(
      And, started_s_ago, solver_->make_term(Not, completed_in_window));
  return solver_->make_term(Not, violated);
}

smt::Term SystemVerilogEncoder::ltl_to_sat(const slang::ast::AssertionExpr & ae,
                                           bool neg,
                                           smt::TermVec & justice)
{
  using namespace slang::ast;

  switch (ae.kind) {
    case AssertionExprKind::Clocking:
      return ltl_to_sat(ae.as<ClockingAssertionExpr>().expr, neg, justice);

    case AssertionExprKind::StrongWeak: {
      auto & sw = ae.as<StrongWeakAssertionExpr>();
      // The strong/weak qualifier only meaningfully differs for a
      // genuine bounded sequence -- must it eventually complete
      // (strong) or not (weak)? Any other shape (already a plain
      // Boolean/temporal expression) is unaffected by the qualifier
      // under this encoder's infinite-lasso semantics; just unwrap.
      if (sw.strength == StrongWeakAssertionExpr::Strong) {
        Term me = match_exists(sw.expr);
        if (me) {
          // strong(seq): a genuine liveness obligation -- the sequence
          // must eventually complete a match.
          return neg ? ltl_make_G(solver_->make_term(Not, me))
                     : ltl_make_F(me, justice);
        }
      }
      return ltl_to_sat(sw.expr, neg, justice);
    }

    case AssertionExprKind::Simple: {
      auto & simple = ae.as<SimpleAssertionExpr>();
      if (auto * named = resolve_named_assertion_ref(simple.expr)) {
        return ltl_to_sat(*named, neg, justice);
      }
      if (simple.repetition) {
        // See the matching check in assertion_expr_to_bool(): route
        // through the general bounded sequence matcher (which throws
        // for an unbounded repeat count) instead of silently ignoring
        // the repetition.
        Term me = match_exists(ae);
        if (!me) return Term();
        return neg ? solver_->make_term(Not, me) : me;
      }
      Term t = expr_to_term(simple.expr);
      if (!t) return Term();
      Term zero = solver_->make_term(0, t->get_sort());
      Term b = solver_->make_term(Distinct, t, zero);
      return neg ? solver_->make_term(Not, b) : b;
    }

    case AssertionExprKind::SequenceConcat: {
      // A bare `##k Q` property has the same truth value as Q under
      // our infinite-time semantics (modulo a front shift).  Unwrap.
      if (auto m = match_const_delay_seq(ae)) {
        return ltl_to_sat(*m->second, neg, justice);
      }
      return Term();
    }

    case AssertionExprKind::Unary: {
      auto & u = ae.as<UnaryAssertionExpr>();
      switch (u.op) {
        case UnaryAssertionOperator::Not:
          return ltl_to_sat(u.expr, !neg, justice);

        case UnaryAssertionOperator::Always:
        case UnaryAssertionOperator::SAlways: {
          // G phi  (positive)  /  !G phi == F !phi  (negated)
          Term phi = ltl_to_sat(u.expr, neg, justice);
          if (!phi) return Term();
          return neg ? ltl_make_F(phi, justice) : ltl_make_G(phi);
        }

        case UnaryAssertionOperator::Eventually:
        case UnaryAssertionOperator::SEventually: {
          // F phi  (positive)  /  !F phi == G !phi  (negated)
          Term phi = ltl_to_sat(u.expr, neg, justice);
          if (!phi) return Term();
          return neg ? ltl_make_G(phi) : ltl_make_F(phi, justice);
        }

        case UnaryAssertionOperator::NextTime:
        case UnaryAssertionOperator::SNextTime: {
          // !X phi == X !phi, so the negation rides along inside phi.
          Term phi = ltl_to_sat(u.expr, neg, justice);
          if (!phi) return Term();
          return ltl_make_X(phi);
        }

        default: return Term();
      }
    }

    case AssertionExprKind::Binary: {
      auto & b = ae.as<BinaryAssertionExpr>();
      switch (b.op) {
        case BinaryAssertionOperator::And:
        case BinaryAssertionOperator::Or: {
          Term l = ltl_to_sat(b.left, neg, justice);
          Term r = ltl_to_sat(b.right, neg, justice);
          if (!l || !r) return Term();
          bool is_and = (b.op == BinaryAssertionOperator::And);
          if (neg) is_and = !is_and;  // De Morgan
          return solver_->make_term(is_and ? And : Or, l, r);
        }

        case BinaryAssertionOperator::Iff: {
          // Children appear in both polarities, so build them
          // positively and apply the outer negation here.  (Temporal
          // operands under `iff` are not negation-normalized; such
          // properties are exotic and out of scope for fairness.)
          Term l = ltl_to_sat(b.left, false, justice);
          Term r = ltl_to_sat(b.right, false, justice);
          if (!l || !r) return Term();
          return neg ? solver_->make_term(Distinct, l, r)
                     : solver_->make_term(Equal, l, r);
        }

        case BinaryAssertionOperator::Implies: {
          // a implies b == !a || b ;  !(a implies b) == a && !b
          Term l = ltl_to_sat(b.left, !neg, justice);
          Term r = ltl_to_sat(b.right, neg, justice);
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
          Term l = ltl_to_sat(b.left, !neg, justice);
          Term r = ltl_to_sat(*rhs, neg, justice);
          if (!l || !r) return Term();
          for (uint32_t i = 0; i < delay; ++i) r = ltl_make_X(r);
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
            Term l = ltl_to_sat(b.left, false, justice);
            Term r = ltl_to_sat(b.right, false, justice);
            if (!l || !r) return Term();
            // until_with: the terminating cycle must also satisfy a.
            Term term = with ? solver_->make_term(And, l, r) : r;
            if (strong) return ltl_make_U(l, term, justice);
            // weak until: a W term == term R (a || term).
            return ltl_make_R(term, solver_->make_term(Or, l, term));
          }
          // Negated, with operands already in negation-normal form
          // (nl = sat(!left), nr = sat(!right)):
          //   !(a U_strong b)      = !a R !b
          //   !(a W b)             = !b U (!a && !b)
          //   !(a U_strong (a&&b)) = !a R (!a || !b)
          //   !(a W (a&&b))        = (!a || !b) U !a
          Term nl = ltl_to_sat(b.left, true, justice);
          Term nr = ltl_to_sat(b.right, true, justice);
          if (!nl || !nr) return Term();
          if (!with) {
            if (strong) return ltl_make_R(nl, nr);
            return ltl_make_U(nr, solver_->make_term(And, nl, nr), justice);
          }
          Term nterm = solver_->make_term(Or, nl, nr);  // !(a && b)
          if (strong) return ltl_make_R(nl, nterm);
          return ltl_make_U(nterm, nl, justice);
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

Term SystemVerilogEncoder::assertion_expr_to_bool(
    const slang::ast::AssertionExpr & ae)
{
  using namespace slang::ast;

  switch (ae.kind) {
    case AssertionExprKind::Clocking: {
      // The clocking event has already been baked into our cycle
      // abstraction; just recurse into the underlying expression.
      return assertion_expr_to_bool(ae.as<ClockingAssertionExpr>().expr);
    }

    case AssertionExprKind::Simple: {
      auto & simple = ae.as<SimpleAssertionExpr>();
      if (auto * named = resolve_named_assertion_ref(simple.expr)) {
        return assertion_expr_to_bool(*named);
      }
      if (simple.repetition) {
        // `expr[*n:m]`/`expr[+]`/`expr[*]`: route through the general
        // bounded sequence matcher (which throws for an unbounded
        // repeat count) instead of silently ignoring the repetition
        // and returning a bare `bool(expr)`.
        return match_exists(ae);
      }
      Term t = expr_to_term(simple.expr);
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
        return assertion_expr_to_bool(*matched->second);
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
      if (Term w = weak_seq_bool(sw.expr)) return w;
      return assertion_expr_to_bool(sw.expr);
    }

    case AssertionExprKind::Unary: {
      auto & u = ae.as<UnaryAssertionExpr>();
      switch (u.op) {
        case UnaryAssertionOperator::Not: {
          Term inner = assertion_expr_to_bool(u.expr);
          if (!inner) return Term();
          return solver_->make_term(Not, inner);
        }
        case UnaryAssertionOperator::Always:
        case UnaryAssertionOperator::SAlways: {
          // Pure safety: `always P` is true at the current cycle
          // exactly when P is true at the current cycle (the
          // "always at every cycle" closure is implicit in the
          // per-cycle property check).
          return assertion_expr_to_bool(u.expr);
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
        Term lhs = assertion_expr_to_bool(*lhs_inner);
        if (!lhs) lhs = match_exists(*lhs_inner);
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
          rhs = assertion_expr_to_bool(*rhs_inner);
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
          if (Term inner = assertion_expr_to_bool(*elem.sequence)) {
            delay += wmax;
            rhs = inner;
            for (uint32_t i = 1; i <= wmax - wmin; ++i) {
              rhs = solver_->make_term(Or, rhs, delay_bool(inner, i));
            }
          }
        } else {
          rhs = assertion_expr_to_bool(*rhs_inner);
        }
        if (!rhs) return Term();

        if (delay > 0) {
          // Materialize the antecedent as a 1-bit BV so the
          // history chain has a value to latch.
          Sort bv1 = solver_->make_sort(BV, 1);
          Term one_bv1 = solver_->make_term(1, bv1);
          Term zero_bv1 = solver_->make_term(0, bv1);
          Term lhs_bv = solver_->make_term(Ite, lhs, one_bv1, zero_bv1);
          Term delayed_bv = make_history_chain(lhs_bv, delay);
          lhs = solver_->make_term(Equal, delayed_bv, one_bv1);
        }
        Term result = solver_->make_term(Implies, lhs, rhs);
        if (lhs_delay > 0) {
          result = solver_->make_term(Or, ltl_before_cycle(lhs_delay), result);
        }
        // `disable iff`: exempt this cycle's check if the disable
        // condition held anywhere in the antecedent-to-consequent
        // shift window, not just at the single cycle the check is
        // anchored at.
        if (Term dw = disable_window(delay)) {
          result = solver_->make_term(Or, dw, result);
        }
        return result;
      }

      Term lhs = assertion_expr_to_bool(b.left);
      Term rhs = assertion_expr_to_bool(b.right);
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

}  // namespace pono
