#include "engines/kinduction.h"
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

// A property expression is evaluated at every clock tick, so
// `assert property (s_eventually q)` means "q infinitely often", not
// "q at least once".  Here q holds at exactly one cycle, so the two
// readings disagree and only the per-cycle one reports the violation.
TEST_P(SVUnitTests, EventuallyNotRecurring)
{
  check_liveness_bmc("eventually_not_recurring.sv", 20);
}

// `disable iff` on a temporal property.  The pair differs only in the
// exemption, so it pins down that the condition is applied at all --
// it used to be computed and then ignored whenever the property did
// not take the safety fast path, making the assertion stronger than
// written.
TEST_P(SVUnitTests, DisableIffTemporalHolds)
{
  check_bmc("disable_iff_temporal.sv", 10, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, DisableIffTemporalFails)
{
  check_bmc("disable_iff_temporal_fails.sv", 1);
}

// The same pairing for a property that stays on the LTL tableau
// path. A liveness attempt's evaluation never ends, so a condition
// that rises after it starts still aborts it; the exemption used to
// cover only the anchor cycle and reported a violation here.
TEST_P(SVUnitTests, DisableIffLivenessHolds)
{
  check_liveness_bmc("disable_iff_liveness.sv", 12, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, DisableIffLivenessFails)
{
  check_liveness_bmc("disable_iff_liveness_fails.sv", 12);
}

// `req |-> ##[1:$] gnt`. An unbounded wait spans no finite window,
// so the property goes to the tableau, where the delay fixes where
// the wait starts and F carries it onward. Each is paired: an F
// whose eventuality is never discharged would satisfy the holds
// case, and one that is vacuously true would satisfy it too.
TEST_P(SVUnitTests, UnboundedDelayHolds)
{
  check_liveness_bmc("unbounded_delay.sv", 12, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, UnboundedDelayFails)
{
  check_liveness_bmc("unbounded_delay_fails.sv", 12);
}

// The minimum has to survive the trip to the tableau: a grant one
// cycle after the request is too early for `##[2:$]`.
TEST_P(SVUnitTests, UnboundedDelayMinimumFails)
{
  check_liveness_bmc("unbounded_delay_min.sv", 12);
}

// `gnt[->n]` and `gnt[=n]`: reaching the n-th occurrence, which is
// an eventuality rather than a finite window. Paired both ways, and
// the count is exercised separately -- a count that collapsed to one
// occurrence would still pass the n=1 pair.
TEST_P(SVUnitTests, GotoRepetitionHolds)
{
  check_liveness_bmc("goto_repetition.sv", 12, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, GotoRepetitionFails)
{
  check_liveness_bmc("goto_repetition_fails.sv", 12);
}

TEST_P(SVUnitTests, GotoRepetitionCountFails)
{
  check_liveness_bmc("goto_repetition_count.sv", 14);
}

TEST_P(SVUnitTests, GotoRepetitionCountHolds)
{
  check_liveness_bmc(
      "goto_repetition_count_holds.sv", 14, ProverResult::UNKNOWN);
}

// A goto or nonconsecutive count as an antecedent. It needs its match
// to end *now*, which no bounded window of offsets can say -- but a
// saturating counter can, and it keeps the implication a plain safety
// property rather than pushing it onto the tableau.

TEST_P(SVUnitTests, GotoRepetitionAntecedent)
{
  check_bmc("goto_repetition_antecedent.sv", 6);
}

TEST_P(SVUnitTests, GotoRepetitionAntecedentEndsOnOccurrence)
{
  check_prover<KInduction>(
      "goto_repetition_antecedent_holds.sv", 20, ProverResult::TRUE);
}

// `[=n]` shares the counter but drops the "ends on an occurrence"
// half, so it fires a cycle later too.
TEST_P(SVUnitTests, NonconsecutiveRepetitionAntecedent)
{
  check_bmc("nonconsec_repetition_antecedent.sv", 7);
}

// weak() over sequence shapes written through different operators.
// Each fixture is the same sequence as the base, so the shared
// refutation depth is the assertion: a misplaced span shifts it or
// drops the refutation.
TEST_P(SVUnitTests, WeakSequenceBase) { check_bmc("weak_seq_fails.sv", 1); }

TEST_P(SVUnitTests, WeakSequenceIntersect)
{
  check_bmc("weak_seq_intersect.sv", 1);
}

TEST_P(SVUnitTests, WeakSequenceThroughout)
{
  check_bmc("weak_seq_throughout.sv", 1);
}

TEST_P(SVUnitTests, WeakSequenceWithin) { check_bmc("weak_seq_within.sv", 1); }

TEST_P(SVUnitTests, WeakSequenceLeadingRepetition)
{
  check_bmc("weak_seq_repetition.sv", 1);
}

// A leading repetition that can match emptily, which leaves no cycle
// marking where an attempt began -- once the only question is which
// attempt is being checked, that no longer matters. Paired with the
// unwrapped spelling at the same depth.
TEST_P(SVUnitTests, WeakSequenceEmptyLeadingRepetition)
{
  check_bmc("weak_seq_empty_leading.sv", 2);
}

TEST_P(SVUnitTests, BareSequenceEmptyLeadingRepetition)
{
  check_bmc("bare_seq_empty_leading.sv", 2);
}

// `and`/`or` over multi-cycle operands, which used to reach the
// tableau and be handed the strong obligation weak withholds. The
// `and` pair is the assertion: operands of spans 1 and 2 sharing a
// start are the 3-cycle chain, so both refute at depth 2.
TEST_P(SVUnitTests, WeakSequenceAnd) { check_bmc("weak_seq_and.sv", 2); }

TEST_P(SVUnitTests, WeakSequenceAndChain)
{
  check_bmc("weak_seq_and_chain.sv", 2);
}

TEST_P(SVUnitTests, WeakSequenceOr) { check_bmc("weak_seq_or.sv", 2); }

// Where `and` places the composite's end, and where it insists the
// operands begin, read off a design whose signals are each true at
// exactly one cycle. The consequent is 1'b0, so the refutation depth
// is the antecedent's own end cycle.
TEST_P(SVUnitTests, SeqAndCommonStart)
{
  check_bmc("seq_and_common_start.sv", 3);
}

// The two negatives, proved rather than left unrefuted: a bounded
// run cannot distinguish "never matches" from "has not matched yet",
// which is exactly what a too-permissive `and` would need it to.
TEST_P(SVUnitTests, SeqAndNoCommonStart)
{
  check_prover<KInduction>(
      "seq_and_no_common_start.sv", 20, ProverResult::TRUE);
}

TEST_P(SVUnitTests, SeqAndOneOperandOnly)
{
  check_prover<KInduction>(
      "seq_and_one_operand_only.sv", 20, ProverResult::TRUE);
}

// `or` over that same rejected pair: a union asks nothing about the
// other operand's start, so it matches where `and` does not.
TEST_P(SVUnitTests, SeqOrIndependentMatches)
{
  check_bmc("seq_or_independent_matches.sv", 3);
}

// ---------------------------------------------------------------------------
// Bounded cycle ranges on the unary property operators
// (`eventually [m:n]`, `s_always [m:n]`, `nexttime [k]`,
// `always [m:$]`).  Each holds/fails pair below differs *only* in the
// window, so the pair pins down that the window is honoured rather
// than dropped -- which it silently was, encoding every one of these
// as if unbounded.
//
// All but the `[m:$]` one are plain safety properties: a bounded
// window names a last cycle, so the check re-anchors there and reads
// backwards.  That is why these use check_bmc(), which additionally
// pins the exact cycle the violation is found at.
// ---------------------------------------------------------------------------

// The regression test for the unsound direction: `count` does reach 3
// repeatedly, so ignoring the [2:3] window leaves the true property
// `s_eventually (count == 3)` and no violation is reported for a
// design that genuinely fails.
TEST_P(SVUnitTests, EventuallyRangeFails)
{
  check_bmc("eventually_range_fails.sv", 3);
}

TEST_P(SVUnitTests, EventuallyRangeHolds)
{
  check_bmc("eventually_range.sv", 12, ProverResult::UNKNOWN);
}

// Nested under an implication, i.e. where the window is relative to
// the match point rather than to the start of the trace.
TEST_P(SVUnitTests, SAlwaysRangeHolds)
{
  check_bmc("s_always_range.sv", 12, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, SAlwaysRangeFails)
{
  check_bmc("s_always_range_fails.sv", 5);
}

// `nexttime [k]` must shift k cycles, not one.  The holds variant used
// to report a spurious counterexample, since a single shift lands on
// count == 1 rather than count == 3.
TEST_P(SVUnitTests, NextTimeRangeHolds)
{
  check_bmc("nexttime_range.sv", 12, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, NextTimeRangeFails)
{
  check_bmc("nexttime_range_fails.sv", 3);
}

// `always [m:$]` -- the windowed form whose upper bound may be
// unbounded, encoded as m forward shifts around the ordinary G tester.
TEST_P(SVUnitTests, AlwaysRangeUnbounded)
{
  check_liveness_bmc("always_range_unbounded.sv", 8, ProverResult::UNKNOWN);
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
      SystemVerilogEncoder::encode(fts, sv_path("multi_assert.sv"));
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
// both set the ConcurrentAssertion handler's `is_cover` flag. A
// multi-cycle sequence goal reaches the same reachability duality as
// CoverProperty, now that a bare sequence is a per-cycle check
// anchored at the attempt's own tick rather than an eventuality.
TEST_P(SVUnitTests, CoverSequence) { check_bmc("cover_sequence.sv", 1); }

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
// ltl_to_sat()'s Simple case recurses into that.
TEST_P(SVUnitTests, NamedSequencePropertyDecl)
{
  check_bmc("named_property_decl.sv", 1);
}

// The same, but parameterized -- how real SVA property libraries are
// written. Each is paired with the property spelled out inline, and
// the shared refutation depth is the assertion: a reference that
// encoded but lost its arguments would not agree with its twin.
TEST_P(SVUnitTests, NamedPropertyArgs)
{
  check_bmc("named_property_args.sv", 1);
}

TEST_P(SVUnitTests, NamedPropertyArgsInline)
{
  check_bmc("named_property_args_inline.sv", 1);
}

// An argument that is a whole subexpression, not a bare signal.
TEST_P(SVUnitTests, NamedPropertyExprArg)
{
  check_bmc("named_property_expr_arg.sv", 1);
}

TEST_P(SVUnitTests, NamedPropertyExprArgInline)
{
  check_bmc("named_property_expr_arg_inline.sv", 1);
}

// A named sequence used directly as a property reaches the tableau
// rather than the safety path, so the substitution has to survive
// there too.
TEST_P(SVUnitTests, NamedSequenceArgs)
{
  check_liveness_bmc("named_sequence_args.sv", 4);
}

TEST_P(SVUnitTests, NamedSequenceArgsInline)
{
  check_bmc("named_sequence_args_inline.sv", 1);
}

// Local variables and recursion still need a binding environment
// that expanding the body does not supply.
TEST_P(SVUnitTests, NamedPropertyLocalVarRejected)
{
  expect_encode_throws("named_property_localvar.sv");
}

// An argument can be a clocking event, which is a route into the
// design's clock that does not look like one syntactically -- the
// multiclock rejection has to see through it.
TEST_P(SVUnitTests, NamedPropertyClockArgRejected)
{
  expect_encode_throws("named_property_clock_arg.sv");
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

// weak(seq) withholds the obligation to match only for the attempts
// still in flight at the end of the trace. Every attempt older than
// the sequence's own maximum span S has had every cycle it could
// ever use, so weak_seq_bool() checks the one that began S cycles
// ago. `a ##1 1'b1` cannot fail after its leading element, which
// makes this exactly as strong as `a`, and `a` is free.
TEST_P(SVUnitTests, WeakSeqLeadingElementRequired)
{
  check_bmc("weak_seq.sv", 1);
}

// The discriminating pair for that: the same sequence and design as
// BareSequenceConcatProperty with `weak` written out, refuted at the
// same depth. Anchoring the check on where the leading element held
// -- rather than on every tick -- made the wrapped form prove while
// the unwrapped one was refuted, though the LRM defines them as the
// same property.
TEST_P(SVUnitTests, WeakSeqMatchesBareSequence)
{
  check_bmc("weak_seq_never_matches.sv", 1);
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
// A named clocking block as an assertion's clocking event. slang
// binds `@(cb)` as a reference to the block rather than to a signal,
// so the clock is the event the block was declared with. (A `default
// clocking` block needs nothing here: slang substitutes its event
// A sampled-value function's `clocking_event` argument naming the
// clock the design already runs on: redundant rather than wrong, so
// it is checked and dropped. Also covers `$past`'s omitted `enable`
// slot, which has to read as "not supplied".
TEST_P(SVUnitTests, SampledValueClockingEvent)
{
  check_prover<KInduction>("sampled_clocking_event.sv", 6, ProverResult::TRUE);
}

// Naming a second clock used to be accepted and ignored, which made
// this design prove that the two clocks agree while a property on
// that same clock was rejected. Both paths now share one record of
// the design's clock.
TEST_P(SVUnitTests, Unsupported_SampledValueSecondClock)
{
  expect_encode_throws("sampled_clocking_event_second_clock.sv");
}

// during elaboration, so the assertion arrives already clocked.)
TEST_P(SVUnitTests, NamedClockingBlockEvent)
{
  check_prover<KInduction>("clocking_block_event.sv", 6, ProverResult::TRUE);
}

// Resolving that reference recurses, so a block on another clock
// reaches the multiclock check rather than slipping past it.
TEST_P(SVUnitTests, Unsupported_ClockingBlockOnSecondClock)
{
  expect_encode_throws("clocking_block_second_clock.sv");
}

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
// Unbounded consecutive sequence repetition. Ending now, a run of at
// least n ends with a run of exactly n, and where the run began
// changes nothing about whether it ends here -- so `[+]` and `[*n:$]`
// come to their bounded twins, which is what each pair below checks
// by sharing a refutation depth.
//
// `[*]` composes like any other repetition, since an empty match is
// reported alongside the offset vector rather than in it. The one
// shape that stays rejected -- a sequence that matches emptily and is
// *all* there is to an antecedent -- is a deliberate non-goal and
// lives in test_systemverilog_unsupported.cpp.
// ---------------------------------------------------------------------------

// A zero lower bound also matches emptily. The empty match absorbs a
// cycle of the following delay, so it is not just "one iteration
// fewer".

TEST_P(SVUnitTests, EmptyRepetitionMatch)
{
  check_prover<KInduction>("empty_repeat_match.sv", 16, ProverResult::TRUE);
}

TEST_P(SVUnitTests, EmptyRepetitionMatchFails)
{
  check_bmc("empty_repeat_match_fails.sv", 2);
}

TEST_P(SVUnitTests, EmptyRepetitionUnbounded)
{
  check_prover<KInduction>("empty_repeat_unbounded.sv", 16, ProverResult::TRUE);
}

// At `##0` an empty match yields no match at all (LRM 16.9.2.1),
// unlike `##1` where it absorbs the delay -- so the empty branch
// drops out rather than acting as an identity.

TEST_P(SVUnitTests, EmptyRepetitionZeroDelay)
{
  check_prover<KInduction>(
      "empty_repeat_zero_delay.sv", 20, ProverResult::TRUE);
}

TEST_P(SVUnitTests, EmptyRepetitionZeroDelayFails)
{
  check_bmc("empty_repeat_zero_delay_fails.sv", 2);
}

// An empty-matching sequence as the *whole* antecedent. It matches at
// every cycle, so the implication is unconditional -- degenerate, but
// well defined, so it encodes with a warning rather than throwing.
// The `Equiv` twin shares its refutation depth with the consequent
// alone, and the `Nonempty` one shows the collapse does not reach a
// repetition that cannot match emptily.

TEST_P(SVUnitTests, EmptyMatchingAntecedent)
{
  check_bmc("unbounded_repeat_star.sv", 6);
}

TEST_P(SVUnitTests, EmptyMatchingAntecedentEquiv)
{
  check_bmc("unbounded_repeat_star_equiv.sv", 6);
}

TEST_P(SVUnitTests, NonemptyRepetitionAntecedentStillGuards)
{
  check_prover<KInduction>(
      "nonempty_repeat_antecedent.sv", 16, ProverResult::TRUE);
}

// Where else an empty-admitting repetition can sit. Each fixture
// states the completion condition worked out by hand and proves the
// encoder never reports a match outside it; the two `Fails` twins
// shift that condition a cycle late, which is the error an empty
// match invites, so none of these can be holding vacuously.

TEST_P(SVUnitTests, EmptyRepetitionTrailing)
{
  check_prover<KInduction>("star_trailing.sv", 18, ProverResult::TRUE);
}

TEST_P(SVUnitTests, EmptyRepetitionTrailingFails)
{
  check_bmc("star_trailing_fails.sv", 4);
}

TEST_P(SVUnitTests, EmptyRepetitionMidSequence)
{
  check_prover<KInduction>("star_middle.sv", 18, ProverResult::TRUE);
}

TEST_P(SVUnitTests, EmptyRepetitionMidSequenceFails)
{
  check_bmc("star_middle_fails.sv", 4);
}

TEST_P(SVUnitTests, EmptyRepetitionDelayRange)
{
  check_prover<KInduction>("star_delay_range.sv", 18, ProverResult::TRUE);
}

TEST_P(SVUnitTests, EmptyRepetitionIntersectOperand)
{
  check_prover<KInduction>("star_intersect.sv", 18, ProverResult::TRUE);
}

TEST_P(SVUnitTests, EmptyRepetitionFirstMatch)
{
  check_prover<KInduction>("star_first_match.sv", 18, ProverResult::TRUE);
}

TEST_P(SVUnitTests, SequenceRepetitionPlus)
{
  check_bmc("unbounded_repeat_plus.sv", 1);
}

TEST_P(SVUnitTests, SequenceRepetitionPlusBounded)
{
  check_bmc("unbounded_repeat_plus_bounded.sv", 1);
}

TEST_P(SVUnitTests, SequenceRepetitionUnboundedRange)
{
  check_bmc("unbounded_repeat_range.sv", 1);
}

TEST_P(SVUnitTests, SequenceRepetitionUnboundedRangeBounded)
{
  check_bmc("unbounded_repeat_range_bounded.sv", 1);
}

// ---------------------------------------------------------------------------
// Property-level connectives ltl_to_sat()/assertion_expr_to_bool()
// previously had no gadget for -- silently dropped (the whole property
// simply never checked, no thrown error) before this session started,
// then converted to a clean throw, and now genuinely implemented by
// composing this file's existing gadgets. `accept_on`/`reject_on` remain
// unimplemented (see Gap_PropertyAcceptOn below) -- their formal
// semantics (the abort condition can supersede the property's outcome
// at *any* cycle during its evaluation, not just at a single recursive
// call) don't localize the way the other four do, so they need their
// own dedicated semantics research before any code is written.
// ---------------------------------------------------------------------------

// In-property `if (sel) a else b`: a plain ITE composition, reduces to
// a current-cycle Boolean via assertion_expr_to_bool()'s own
// Conditional case, so this is an ordinary safety property (not routed
// through the LTL/justice machinery) -- sel/a/b are all free, so BMC
// finds a violation (e.g. sel=1, a=0) immediately.
TEST_P(SVUnitTests, PropertyConditional)
{
  check_bmc("property_conditional.sv", 0);
}

// In-property `case (sel) 0: a; 1: b; default: 1'b1; endcase`: same
// current-cycle-Boolean fast path as PropertyConditional, generalized
// to N branches -- sel/a/b free, BMC finds a violation (e.g. sel=0,
// a=0) immediately.
TEST_P(SVUnitTests, PropertyCase) { check_bmc("property_case.sv", 0); }

// `accept_on`/`reject_on`/`sync_accept_on`/`sync_reject_on`: deferred,
// see the file-level comment above this section.
TEST_P(SVUnitTests, Gap_PropertyAcceptOn)
{
  check_bmc("property_accept_on.sv", 1, ProverResult::UNKNOWN);
}

// `a intersect b` used directly as a property (as opposed to as the
// antecedent of `|->`/`|=>`, where offsets_ending_now() already
// handles it -- see SeqIntersect above): per the LRM a bare sequence
// used as a property has implicit `strong` semantics ("must eventually
// match"), so this routes through try_strong_sequence()'s match_exists()
// + make_F() -- a genuine liveness obligation, hence check_liveness_bmc()
// rather than check_bmc(). `a`/`b` are both free every cycle, so a
// trace where they're never simultaneously true (e.g. always a=0)
// violates the obligation.
TEST_P(SVUnitTests, PropertyIntersectTopLevel)
{
  check_bmc("property_intersect_toplevel.sv", 0);
}

// `a #-# b`: a required (not merely conditional) sequential
// composition -- match_exists(a) (trivially "a holds now" for a
// length-1 sequence) AND b, checked at every cycle. Reduces to a plain
// current-cycle Boolean in principle, but this encoder doesn't (yet)
// give FollowedBy an assertion_expr_to_bool() fast path the way
// Conditional/Case got, so it's still routed through the LTL/justice
// machinery -- check_liveness_bmc() rather than check_bmc(). `a`/`b`
// free, so BMC finds a=0 (or b=0) immediately.
TEST_P(SVUnitTests, PropertyFollowedBy)
{
  check_liveness_bmc("property_followed_by.sv", 3);
}

// A bare multi-element sequence (`a ##1 b`) used directly as a
// property: same implicit-`strong` reasoning as
// PropertyIntersectTopLevel -- must eventually match. `a`/`b` free, so
// a trace that never has a followed by b one cycle later violates it.
TEST_P(SVUnitTests, BareSequenceConcatProperty)
{
  check_bmc("bare_sequence_concat_property.sv", 1);
}

// Temporal (non-safety) `assume`/`restrict property`: a fairness
// constraint. The assumption's own tableau is required at every
// cycle, and the eventualities it may promise and never keep become
// justice conditions on the counterexample lasso, alongside the
// property's own.
//
// Each of these is proved rather than merely left unrefuted. BMC can
// only fail to find a lasso, which looks identical whether the
// assumption bites or was dropped on the floor.
TEST_P(SVUnitTests, TemporalAssumeProperty)
{
  check_liveness_prover<KInduction>(
      "temporal_assume_property.sv", 20, ProverResult::TRUE);
}

TEST_P(SVUnitTests, TemporalAssumeAbsent)
{
  check_liveness_bmc("temporal_assume_absent.sv", 5);
}

// The negative that keeps the positives honest: an assumption whose
// tableau could not be satisfied would empty the model and prove
// everything, including this.
TEST_P(SVUnitTests, TemporalAssumeUnrelated)
{
  check_liveness_bmc("temporal_assume_unrelated.sv", 5);
}

TEST_P(SVUnitTests, TemporalRestrictProperty)
{
  check_liveness_prover<KInduction>(
      "temporal_restrict_property.sv", 20, ProverResult::TRUE);
}

// The construct's actual use: an acknowledgement the design cannot
// force and the environment has to promise.
TEST_P(SVUnitTests, FairnessHandshake)
{
  check_liveness_prover<KInduction>(
      "fairness_handshake.sv", 20, ProverResult::TRUE);
}

TEST_P(SVUnitTests, FairnessHandshakeUnfair)
{
  check_liveness_bmc("fairness_handshake_unfair.sv", 5);
}

// An assumption constrains the assertions written before it too,
// which is why the conditions are distributed after the walk.
TEST_P(SVUnitTests, FairnessAssumeAfterAssert)
{
  check_liveness_prover<KInduction>(
      "fairness_assume_after_assert.sv", 20, ProverResult::TRUE);
}

TEST_P(SVUnitTests, FairnessTwoAssumptions)
{
  check_liveness_prover<KInduction>(
      "fairness_two_assumptions.sv", 20, ProverResult::TRUE);
}

// `restrict` is an assumption too, and shares the whole path with
// `assume`; here both kinds pool their conditions, and the assertion
// is shaped so only the `restrict`'s can carry it.
TEST_P(SVUnitTests, FairnessAssumeAndRestrict)
{
  check_liveness_prover<KInduction>(
      "fairness_assume_and_restrict.sv", 20, ProverResult::TRUE);
}

// Where a fairness constraint stops, shown on one design asked both
// ways. No finite trace contradicts `s_eventually a`, so the safety
// property is refuted; the liveness one, which is about infinite
// traces, is proved. Same design, same assumption -- and the split
// BTOR2 already makes between `fair` and `bad`.
TEST_P(SVUnitTests, FairnessNotAppliedToSafety)
{
  check_bmc("fairness_not_applied_to_safety.sv", 1);
}

TEST_P(SVUnitTests, FairnessAppliedToLiveness)
{
  check_liveness_prover<KInduction>(
      "fairness_applied_to_liveness.sv", 20, ProverResult::TRUE);
}

TEST_P(SVUnitTests, Unsupported_TemporalAssumeDisableIff)
{
  expect_encode_throws("temporal_assume_disable_iff.sv");
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVSvaTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

// ---------------------------------------------------------------------------
// Consecutive repetition of a whole sequence, `(seq)[*n]`, which the
// LRM makes `seq` concatenated with itself at `##1` -- so its offsets
// are `seq`'s own convolved with themselves.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, SequenceRepetition)
{
  check_prover<KInduction>("sequence_repetition.sv", 20, ProverResult::TRUE);
}

TEST_P(SVUnitTests, SequenceRepetitionFails)
{
  check_bmc("sequence_repetition_fails.sv", 5);
}

TEST_P(SVUnitTests, SequenceRepetitionRange)
{
  check_prover<KInduction>(
      "sequence_repetition_range.sv", 20, ProverResult::TRUE);
}

// ---------------------------------------------------------------------------
// Sequences whose match starts an unbounded distance back -- an
// unbounded inter-element delay, or a counted repetition composed
// with another element. No window of offsets spans one, but where it
// *ends* is still a definite cycle, which is all an implication's
// antecedent needs. Elsewhere such a sequence is an eventuality, and
// the matcher declines so the tableau models it instead.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, UnboundedDelayAntecedent)
{
  check_prover<KInduction>(
      "unbounded_delay_antecedent.sv", 20, ProverResult::TRUE);
}

TEST_P(SVUnitTests, UnboundedDelayAntecedentFails)
{
  check_bmc("unbounded_delay_antecedent_fails.sv", 3);
}

TEST_P(SVUnitTests, GotoThenElementAntecedent)
{
  check_prover<KInduction>("goto_then_antecedent.sv", 20, ProverResult::TRUE);
}

// In a consequent or as a bare property, an unbounded-span sequence
// is an eventuality rather than a safety obligation: the matcher
// declines and the tableau carries it, with "a match ends at this
// cycle" as the condition F must discharge.

TEST_P(SVUnitTests, UnboundedDelayConsequentFails)
{
  check_liveness_bmc("unbounded_delay_consequent.sv", 4);
}

TEST_P(SVUnitTests, UnboundedDelayConsequentHolds)
{
  check_liveness_bmc(
      "unbounded_delay_consequent_holds.sv", 12, ProverResult::UNKNOWN);
}

// A counted repetition *following* another element counts from where
// that element ended, not from the beginning of time -- one register
// per occurrence count, shifted as the count rises.

TEST_P(SVUnitTests, GotoAfterElement)
{
  check_prover<KInduction>("goto_after_element.sv", 40, ProverResult::TRUE);
}

TEST_P(SVUnitTests, GotoAfterElementFails)
{
  check_bmc("goto_after_element_fails.sv", 9);
}

TEST_P(SVUnitTests, GotoAfterElementDelay)
{
  check_prover<KInduction>(
      "goto_after_element_delay.sv", 40, ProverResult::TRUE);
}

// At `##0` the window opens on the prefix's own last cycle, so a
// record made there belongs to the count before that cycle's
// occurrence -- and the window may be the one opening now, which no
// register has seen yet.

TEST_P(SVUnitTests, GotoAfterElementOverlap)
{
  check_prover<KInduction>(
      "goto_after_element_overlap.sv", 40, ProverResult::TRUE);
}

TEST_P(SVUnitTests, GotoAfterElementOverlapFails)
{
  check_bmc("goto_after_element_overlap_fails.sv", 4);
}

}  // namespace pono_tests
