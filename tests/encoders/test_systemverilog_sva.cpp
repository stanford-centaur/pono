#ifdef WITH_SLANG

#include "sv_test_fixture.h"

using namespace pono;
using namespace smt;

namespace pono_tests {

// ---------------------------------------------------------------------------
// Migrated from the original test_systemverilog.cpp -- already compositional
// and (for the ReqAck pair) already the holds+fails template the rubric
// asks every property-style test to generalize; kept as-is per the triage
// in the plan.
// ---------------------------------------------------------------------------

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
  SystemVerilogEncoder enc(sv_path("multi_assert.sv"), fts);
  EXPECT_EQ(enc.propvec().size(), 3u);
}

// ---------------------------------------------------------------------------
// TOP PRIORITY (see the plan's Context section): `assume`/`restrict
// property` are a complete silent no-op today -- not even a log line --
// which is the single most dangerous gap this suite documents, since a
// dropped environment assumption can make a real bug look like a proof.
// `cover property` dropping silently is comparatively low-stakes (no proof
// obligation to corrupt), so it's checked for clean absence instead of
// being asserted against a safety property.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, Gap_AssumePropertyConstrainsTrace)
{
  check_bmc("assume_property.sv", 4, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, Gap_RestrictPropertyConstrainsTrace)
{
  check_bmc("restrict_property.sv", 4, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, CoverPropertyCleanlyAbsent)
{
  expect_silently_dropped("cover_property.sv");
}

// ---------------------------------------------------------------------------
// $rose/$fell/$changed/$onehot/$onehot0/$isunknown, checked as identities.
// ---------------------------------------------------------------------------

// GAP (confirmed empirically): only $past/$stable are handled in
// expr_to_term()'s Call case; $rose throws "unsupported call to
// $rose" before $fell/$changed are even reached.
TEST_P(SVUnitTests, Gap_RoseFellChangedHold)
{
  check_bmc("rose_fell_changed.sv", 4, ProverResult::UNKNOWN);
}

// GAP (confirmed empirically, after fixing an unrelated syntax bug in
// an earlier draft that mixed a property-level `|->` inside a boolean
// `&&` expression): $onehot throws "unsupported call to $onehot"
// before $onehot0/$isunknown are even reached -- the same "only
// $past/$stable are handled" gap as Gap_RoseFellChangedHold above.
TEST_P(SVUnitTests, Gap_OnehotIsUnknownHold)
{
  check_bmc("onehot_isunknown.sv", 4, ProverResult::UNKNOWN);
}

// ---------------------------------------------------------------------------
// Range delay (##[m:n]), first_match, and a named sequence/property
// declaration referenced by name.  GAP (confirmed empirically): all three
// leave propvec()/ltl_justice() completely empty -- the whole assertion is
// silently dropped, the same dangerous contract as assume/cover/restrict
// above, just reached through a different unsupported construct.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, Gap_RangeDelaySeq)
{
  expect_silently_dropped("range_delay.sv");
}

TEST_P(SVUnitTests, Gap_FirstMatchSeq)
{
  expect_silently_dropped("first_match_seq.sv");
}

// GAP (confirmed empirically): referencing a named `property`
// declaration throws "unsupported expression kind 39" -- a
// different failure mode again (an uncaught exception, not a silent
// drop), since resolving the AssertionInstanceExpression that names
// p_check apparently routes through expr_to_term() rather than
// assertion_expr_to_bool()/ltl_to_sat().
TEST_P(SVUnitTests, Gap_NamedSequencePropertyDecl)
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
// GAP (confirmed empirically): strong(seq)/weak(seq) applied to a plain
// SequenceConcat leave both propvec() and ltl_justice() empty -- despite
// AssertionExprKind::StrongWeak appearing in ltl_to_sat()'s switch, this
// particular combination falls through silently rather than producing a
// liveness obligation (strong) or a vacuously-true safety property (weak).
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, Gap_StrongSeqObligation)
{
  expect_silently_dropped("strong_seq.sv");
}

TEST_P(SVUnitTests, Gap_WeakSeqVacuousHold)
{
  expect_silently_dropped("weak_seq.sv");
}

// ---------------------------------------------------------------------------
// Sequence intersect/within/throughout: absent from the AssertionExprKind
// switches (confirmed by inspection), so the whole assertion is silently
// dropped rather than thrown.  These are lower-severity than the
// assume/cover/range-delay/strong/weak cases above only in the sense that
// they were never expected to be in scope yet; the silent-drop *mechanism*
// is identical.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, SeqIntersectSilentlyDropped)
{
  expect_silently_dropped("seq_intersect.sv");
}

TEST_P(SVUnitTests, SeqWithinSilentlyDropped)
{
  expect_silently_dropped("seq_within.sv");
}

TEST_P(SVUnitTests, SeqThroughoutSilentlyDropped)
{
  expect_silently_dropped("seq_throughout.sv");
}

// ---------------------------------------------------------------------------
// GAP (confirmed empirically): a multiclock property is not "honored" or
// "cleanly rejected" or even "mis-encoded as single-clock" -- it is
// silently dropped entirely (propvec()/ltl_justice() both empty), most
// likely because the nested mid-sequence clock change doesn't match the
// shape match_const_delay_seq()/the Clocking-unwrap loop expect.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, Gap_MulticlockProperty)
{
  expect_silently_dropped("multiclock_property.sv");
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVSvaTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests

#endif
