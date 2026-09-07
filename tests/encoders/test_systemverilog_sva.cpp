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

// `iff` with a temporal (eventuality) operand on one side, forcing
// the whole property through ltl_to_sat()'s Iff case -- previously
// that case always visited its operands positively regardless of the
// requested polarity, so a temporal operand never got its correct
// negation-normalized dual construction. `done` is tied low and
// `req` tied high, so `(s_eventually done) iff req` is a constant
// false -- violated immediately.
TEST_P(SVUnitTests, IffTemporal) { check_liveness_bmc("iff_temporal.sv", 2); }

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

// Plain-boolean `iff`, exercising assertion_expr_to_bool()'s safety
// fast path (not ltl_to_sat() -- see IffTemporal below for that).
TEST_P(SVUnitTests, BinaryIff) { check_bmc("binary_iff.sv", 0); }

TEST_P(SVUnitTests, PastCall) { check_bmc("past_call.sv", 1); }

// `$past(expr, n, enable)`'s `enable` argument gates whether each cycle's
// sample is taken -- freezing the whole delayed history on disabled
// cycles rather than being silently ignored (as it was before this
// history-chain gating was added). Checked as a genuine identity against
// a hand-rolled "sample-on-enable, else hold" register that implements
// the same rule directly -- see the comment in past_call_enable.sv.
TEST_P(SVUnitTests, PastCallWithEnable)
{
  check_bmc("past_call_enable.sv", 4, ProverResult::UNKNOWN);
}

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

// `cover sequence(S)` is treated the same as `cover property(P)` --
// both set the ConcurrentAssertion handler's `is_cover` flag. Since
// `a ##1 b` is a genuinely multi-cycle sequence, it hits the
// temporal/sequence-shaped cover-goal throw (same reachability-duality
// contract as CoverProperty above, once implemented): extending
// reachability duality through the LTL tableau for cover goals is a
// real gap, not a deliberate non-goal.
TEST_P(SVUnitTests, Gap_CoverSequence) { check_bmc("cover_sequence.sv", 1); }

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
// 2-element SequenceConcat with a mid-sequence clock change. This
// encoder has no clock-domain-crossing model (no clock dividers, no
// nondeterministic per-cycle choice of which clock toggles), so
// rather than silently collapsing `clk2` onto the same global cycle
// as `clk1` (the property's own outer clock, established first),
// check_clock() rejects the design outright -- correctly out of
// scope, not a gap to eventually close.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, MulticlockPropertyRejected)
{
  expect_encode_throws("multiclock_property.sv");
}

// Same signal, opposite edges (`@(posedge clk)` on one property,
// `@(negedge clk)` on another): check_clock() tracks (signal, edge)
// pairs, not just the signal, so this is rejected for the same reason
// as MulticlockPropertyRejected above.
TEST_P(SVUnitTests, MixedEdgePropertyRejected)
{
  expect_encode_throws("mixed_edge_property.sv");
}

// ---------------------------------------------------------------------------
// Unbounded consecutive sequence repetition (`[*]`, `[+]`, `[*n:$]`) is
// mainstream verification-relevant SVA (e.g. `req[*1:$] ##1 gnt`
// handshake idioms), not a deliberate non-goal -- unlike a truly
// unbounded `forever` loop, an unbounded repeat has an obvious bounded
// approximation (unroll up to the BMC bound), so this is a real,
// worth-fixing gap in offsets_ending_now()'s compile-time-bounded
// model, not an inherent modeling impossibility. It currently throws
// a clear error rather than silently dropping the assertion; encoding
// any partial/approximate result to check a property against would
// require implementing that bounded-unrolling approximation first, so
// these stay throw-based for now (see Gap_CoverSequence in
// test_systemverilog_unsupported.cpp for the same tradeoff).
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, Gap_SequenceRepetitionStar)
{
  expect_encode_throws("unbounded_repeat_star.sv");
}

TEST_P(SVUnitTests, Gap_SequenceRepetitionPlus)
{
  expect_encode_throws("unbounded_repeat_plus.sv");
}

TEST_P(SVUnitTests, Gap_SequenceRepetitionUnboundedRange)
{
  expect_encode_throws("unbounded_repeat_range.sv");
}

// ---------------------------------------------------------------------------
// Property-level connectives ltl_to_sat()/assertion_expr_to_bool() have no
// gadget for -- previously silently dropped (the whole property simply
// never checked, no thrown error) rather than throwing; now throws a
// clear error naming the unsupported shape. Each is mainstream
// verification-relevant SVA, not a deliberate non-goal: in-property
// if/case could plausibly ITE-compose already-built Booleans;
// accept_on/reject_on could plausibly reuse the disable_window()
// machinery `disable iff` already has; intersect/within/throughout/
// followed-by as top-level connectives (as opposed to inside a bounded
// sequence match, which offsets_ending_now() already handles -- see
// SeqIntersect/SeqWithin/SeqThroughout above) and a bare multi-element
// sequence used directly as a property could plausibly delegate to
// offsets_ending_now()/match_exists() and the existing F/G/U/R tableau
// gadgets. None of that plausible follow-up work is attempted here, so
// (matching Gap_UserFunctionCall's convention) these assert UNKNOWN as
// a placeholder for "the correct answer, once implemented" rather than
// a specifically-reasoned verdict.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, Gap_PropertyConditional)
{
  check_bmc("property_conditional.sv", 1, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, Gap_PropertyCase)
{
  check_bmc("property_case.sv", 1, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, Gap_PropertyAcceptOn)
{
  check_bmc("property_accept_on.sv", 1, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, Gap_PropertyIntersectTopLevel)
{
  check_bmc("property_intersect_toplevel.sv", 1, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, Gap_PropertyFollowedBy)
{
  check_bmc("property_followed_by.sv", 1, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, Gap_BareSequenceConcatProperty)
{
  check_bmc("bare_sequence_concat_property.sv", 1, ProverResult::UNKNOWN);
}

// Temporal (non-safety) `assume`/`restrict property` -- previously
// silently dropped (logged and skipped, with the model left less
// constrained than the source describes, risking a spurious
// counterexample from a later `assert`); now throws instead. A real
// gap (the dual of the justice-based proving machinery already built
// for `assert`), not an inherent impossibility.
TEST_P(SVUnitTests, Gap_TemporalAssumeProperty)
{
  check_bmc("temporal_assume_property.sv", 1, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, Gap_TemporalRestrictProperty)
{
  check_bmc("temporal_restrict_property.sv", 1, ProverResult::UNKNOWN);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVSvaTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
