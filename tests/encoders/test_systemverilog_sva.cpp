#ifdef WITH_SLANG

#include "sv_test_fixture.h"

using namespace pono;
using namespace smt;

namespace pono_tests {

TEST_P(SVUnitTests, ReqAckLiveness) { check_liveness_bmc("req_ack.sv", 10); }

TEST_P(SVUnitTests, ReqAckHolds)
{
  check_liveness_bmc("req_ack_holds.sv", 8, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, UntilLiveness) { check_liveness_bmc("until_live.sv", 10); }

TEST_P(SVUnitTests, EventuallyAssertion)
{
  check_liveness_bmc("eventually_assertion.sv", 5);
}

TEST_P(SVUnitTests, AlwaysAssertion) { check_bmc("always_assertion.sv", 5); }

TEST_P(SVUnitTests, BinaryImplication)
{
  check_bmc("binary_implication.sv", 1);
}

TEST_P(SVUnitTests, BinaryNonOverlap) { check_bmc("binary_nonoverlap.sv", 1); }

TEST_P(SVUnitTests, BinaryAnd) { check_bmc("binary_and.sv", 0); }

TEST_P(SVUnitTests, PastCall) { check_bmc("past_call.sv", 1); }

TEST_P(SVUnitTests, SequenceDelay) { check_bmc("sequence_delay.sv", 2); }

TEST_P(SVUnitTests, MultipleAssertions)
{
  SmtSolver s = create_solver(GetParam());
  FunctionalTransitionSystem fts(s);
  auto sv_result =
      SystemVerilogEncoder::encode(sv_path("multi_assert.sv"), fts);
  EXPECT_EQ(sv_result.propvec.size(), 3u);
}

// ---------------------------------------------------------------------------
// `assume property (P)`/`restrict property (P)`: the ConcurrentAssertion
// handler routes the same property-shape computation used for `assert`
// through fts_.add_constraint() instead of propvec_, for the safety
// (non-temporal) fast path -- so an assumed property actually constrains
// every reachable trace rather than being ignored. `cover property` is
// covered separately -- see CoverProperty below.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, AssumePropertyConstrainsTrace)
{
  check_bmc("assume_property.sv", 4, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, RestrictPropertyConstrainsTrace)
{
  check_bmc("restrict_property.sv", 4, ProverResult::UNKNOWN);
}

// `cover property (P)`: implemented via reachability duality (see the
// design-decision note at the top of
// frontends/systemverilog/assertion_walker.cpp) --
// checked exactly like `assert property (!P)`, so finding a
// counterexample to that surrogate assertion is precisely "P was
// reached". `data` is free, so the cover point (data == 5) is reachable
// at the earliest possible cycle.
TEST_P(SVUnitTests, CoverProperty) { check_bmc("cover_property.sv", 1); }

// Procedural immediate `cover(expr);` (distinct from concurrent `cover
// property (...)` above), same reachability-duality contract.
TEST_P(SVUnitTests, ImmediateCover) { check_bmc("immediate_cover.sv", 1); }

// ---------------------------------------------------------------------------
// $rose/$fell/$changed/$onehot/$onehot0/$isunknown. $rose/$fell/$changed
// each build their own 1-cycle latch chain via the same make_history_chain()
// helper $past/$stable use; $onehot/$onehot0 are the standard
// (x & (x-1)) == 0 power-of-two bit trick; $isunknown is always false,
// since this encoder's pure 2-valued bitvector model has no X/Z
// representation at all. Checked as identities here, not one hand-picked
// value.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, RoseFellChangedHold)
{
  check_bmc("rose_fell_changed.sv", 4, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, OnehotIsUnknownHold)
{
  check_bmc("onehot_isunknown.sv", 4, ProverResult::UNKNOWN);
}

// ---------------------------------------------------------------------------
// Range delay (`##[m:n]`) on a single-element consequent: generalized
// into an OR over "consequent held i cycles ago" for i spanning the
// window, anchored at the window's *latest* cycle (the earliest point a
// violation becomes certain). Semantics: whenever `arm` holds at cycle
// i, `data == 10` must hold at *some* cycle in [i+1, i+3]. `arm`/`data`
// are free from cycle 0 (neither is reset-gated), so BMC falsifies by
// holding arm at cycle 0 and data != 10 at cycles 1, 2, and 3 -- one
// cycle earlier than a naive "everything starts at cycle 1" guess would
// suggest, since only registered/output state gets that one-cycle reset
// delay.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, RangeDelaySeq) { check_bmc("range_delay.sv", 3); }

// `first_match(seq)` restricts a (possibly multi-match) sequence to its
// earliest match, which doesn't change whether the sequence matches at
// all -- unwrapped by offsets_ending_now() (first_match doesn't change
// match existence, only which match is *reported*, which this encoder
// never needs to distinguish). `first_match(a ##[1:2] b) |-> 1'b0` is
// falsified as soon as the sequence matches at all; earliest match is
// a@0, b@1 (both free from cycle 0), violated at cycle 1.
TEST_P(SVUnitTests, FirstMatchSeq) { check_bmc("first_match_seq.sv", 1); }

// A named `sequence`/`property` declaration referenced by name:
// `assert property (p_check);` binds p_check as a SimpleAssertionExpr
// wrapping an AssertionInstanceExpression, not a plain boolean
// Expression. Slang has already expanded the referenced property's own
// body (clocking, |->, etc. intact) into
// AssertionInstanceExpression::body; assertion_expr_to_bool()/
// ltl_to_sat()'s Simple case recurses into that, scoped to the
// no-argument, non-recursive case (an argumented or recursive
// property/sequence instantiation throws a clear error instead).
TEST_P(SVUnitTests, NamedSequencePropertyDecl)
{
  check_bmc("named_property_decl.sv", 1);
}

// ---------------------------------------------------------------------------
// The genuinely-holding side of an until-family property, exercising the
// Release tableau gadget (ltl_make_R) that SVA has no direct keyword for.
// Passes.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, UntilHoldsViaReleaseTableau)
{
  check_liveness_bmc("until_holds.sv", 6, ProverResult::UNKNOWN);
}

// ---------------------------------------------------------------------------
// strong(seq)/weak(seq) applied to a general (multi-element) sequence --
// `a ##1 1'b1` is a 2-element SequenceConcat, handled by the general
// sequence matcher (offsets_ending_now()), not just the single-element
// match_const_delay_seq() fast path.
//
// strong(seq) is a genuine liveness obligation: the sequence must
// *eventually* match. Built as `F(match_exists(seq))` in ltl_to_sat()'s
// StrongWeak case. `a` is free and can stay low forever, so BMC (via
// the L2S translator) finds a fair lasso where the eventuality is never
// discharged.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, StrongSeqObligation)
{
  check_liveness_bmc("strong_seq.sv", 5);
}

// weak(seq) carries no obligation to ever match, but an attempt that
// *did* begin must not be a definite, provable failure -- see
// weak_seq_bool(): an attempt began exactly S cycles ago (S = the
// sequence's own max span) and nothing completed anywhere in that
// window. Here, `a` never firing means no attempt ever began at all
// (vacuously fine), and once it does fire, the continuation `##1 1'b1`
// is an unconditional truth that can never itself fail -- so the
// property holds vacuously forever, a genuine tautology (confirmed by
// weak_seq_bool()'s formula reducing to a logical contradiction for
// this shape), not merely "unproven within this bound".
TEST_P(SVUnitTests, WeakSeqVacuousHold)
{
  check_bmc("weak_seq.sv", 4, ProverResult::UNKNOWN);
}

// ---------------------------------------------------------------------------
// Sequence intersect/within/throughout, built directly on
// offsets_ending_now()'s offset vector: Intersect ANDs the two
// operands' vectors entry-by-entry (same span required); Within ORs
// the antecedent's merged "matches here" term over the consequent's own
// window (window_or()); Throughout ANDs a plain boolean over that same
// window (window_and()). Each fixture wraps its sequence composition in
// `|-> 1'b0` (the idiom first_match_seq.sv uses) so violation happens
// as soon as the composite sequence matches at all; see each .sv file
// for the match semantics. `a`/`b`/`c` are free from cycle 0, not just
// from cycle 1, which is why each violates one cycle earlier than a
// naive "everything starts at cycle 1" guess would suggest.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, SeqIntersect) { check_bmc("seq_intersect.sv", 1); }

TEST_P(SVUnitTests, SeqWithin) { check_bmc("seq_within.sv", 1); }

TEST_P(SVUnitTests, SeqThroughout) { check_bmc("seq_throughout.sv", 1); }

// ---------------------------------------------------------------------------
// A multiclock property: the antecedent `a ##1 @(posedge clk2) b` is a
// 2-element SequenceConcat with a mid-sequence clock change.
// offsets_ending_now() unwraps a nested Clocking node exactly like the
// outer handler already does for the top-level one -- per this file's
// documented multiclock design decision, every named clock is treated
// as the same global pono-cycle. This reduces the property to exactly
// the same shape as a plain `a ##1 b |-> 1'b0`. `a`/`b` are free from
// cycle 0 (neither is reset-gated); earliest match is a@0, b@1,
// violated at cycle 1. A real clock-domain-crossing model (relative
// clock ratios, distinct per-domain sampling) would be a substantially
// larger feature; this is the minimal correct behavior given the
// encoder's existing single-clock assumption.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, MulticlockProperty)
{
  check_bmc("multiclock_property.sv", 1);
}

// ---------------------------------------------------------------------------
// Unbounded consecutive sequence repetition (`[*]`, `[+]`, `[*n:$]`): a
// genuine architectural boundary of offsets_ending_now()'s compile-
// time-bounded model, not a "not implemented yet" gap -- an unbounded
// repeat count can't be unrolled into a finite offset vector. Must
// throw a clear error rather than silently dropping the assertion (the
// same contract as runtime-dependent break/continue/while conditions
// elsewhere in this encoder).
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, Unsupported_SequenceRepetitionStar)
{
  expect_encode_throws("unbounded_repeat_star.sv");
}

TEST_P(SVUnitTests, Unsupported_SequenceRepetitionPlus)
{
  expect_encode_throws("unbounded_repeat_plus.sv");
}

TEST_P(SVUnitTests, Unsupported_SequenceRepetitionUnboundedRange)
{
  expect_encode_throws("unbounded_repeat_range.sv");
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVSvaTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests

#endif
