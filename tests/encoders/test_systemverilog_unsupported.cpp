#ifdef WITH_SLANG

#include "sv_test_fixture.h"

using namespace pono;
using namespace smt;

namespace pono_tests {

// ---------------------------------------------------------------------------
// Ledger of genuinely out-of-scope IEEE 1800-2017 constructs (OOP/classes,
// randomization, DPI, functional coverage, programs/checkers/specify,
// fork/join/wait/force-release, non-integral types, dynamic containers,
// defparam/bind). Each is checked via either expect_encode_throws() or
// expect_encode_succeeds_ignoring(), whichever matches how the encoder
// actually rejects it; see the per-test comment when that isn't obvious
// from the test name. Tests prefixed `Gap_` instead cover mainstream-RTL
// features or lvalue-resolution edge cases this encoder doesn't yet
// support -- missed synthesizable-subset work, not deliberate non-goals.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, Unsupported_ClassDecl)
{
  expect_encode_throws("class_decl.sv");
}

TEST_P(SVUnitTests, Unsupported_RandomizeConstraint)
{
  expect_encode_throws("randomize_constraint.sv");
}

TEST_P(SVUnitTests, Unsupported_DpiImport)
{
  expect_encode_throws("dpi_import.sv");
}

TEST_P(SVUnitTests, Unsupported_CovergroupDecl)
{
  expect_encode_throws("covergroup_decl.sv");
}

// `program`/`checker` instances are a verification-only construct with
// no functional-logic counterpart: process_instance() recognizes a
// program instance via DefinitionKind::Program and skips it, and a
// checker instance is a distinct SymbolKind::CheckerInstance the usual
// member walk doesn't match at all -- either way logged via
// logger.log(1, "... ignoring ... instance ...") rather than thrown,
// per the "simulation-only constructs are dropped and logged" half of
// the constructor's documented contract (see
// SystemVerilogEncoder::SystemVerilogEncoder()'s doc comment).
TEST_P(SVUnitTests, Unsupported_ProgramBlock)
{
  expect_encode_succeeds_ignoring("program_block.sv");
}

TEST_P(SVUnitTests, Unsupported_CheckerBlock)
{
  expect_encode_succeeds_ignoring("checker_block.sv");
}

// `fork`/`join` and `wait` are simulation-timing constructs with no
// per-cycle counterpart in this encoder's model; process_statement()'s
// default case logs a warning (logger.log(1, "... skipping unsupported
// statement kind ...")) and skips them, rather than throwing.
TEST_P(SVUnitTests, Unsupported_ForkJoin)
{
  expect_encode_succeeds_ignoring("fork_join.sv");
}

// `case (x) matches ... endcase` (StatementKind::PatternCase) is a
// distinct statement kind from plain case/casex/casez
// (StatementKind::Case) that pre_scan_state_vars()'s
// collect_blocking_targets()/collect_nonblocking_targets() don't
// recognize either -- but since process_statement() itself also
// doesn't process it (falling to the generic unhandled-statement-kind
// default below), the two omissions are consistent: no write inside
// it is ever pre-scanned *or* applied, and the skip is logged.
TEST_P(SVUnitTests, Gap_PatternCase)
{
  expect_encode_succeeds_ignoring("pattern_case.sv");
}

TEST_P(SVUnitTests, Unsupported_WaitStmt)
{
  expect_encode_succeeds_ignoring("wait_stmt.sv");
}

// `expect (property_expr);` is a procedural blocking-wait statement
// (pause until the property holds), not a checked invariant -- the
// same simulation-only category as `wait` above.
TEST_P(SVUnitTests, Unsupported_ExpectProperty)
{
  expect_encode_succeeds_ignoring("expect_property.sv");
}

// `cover sequence(S)` is treated the same as `cover property(P)` --
// both set the ConcurrentAssertion handler's `is_cover` flag. Since
// `a ##1 b` is a genuinely multi-cycle sequence, it hits the
// temporal/sequence-shaped cover-goal throw (same as
// `cover property (a ##1 b)` already does): extending reachability
// duality through the LTL tableau for cover goals is out of scope.
TEST_P(SVUnitTests, Gap_CoverSequence)
{
  expect_encode_throws("cover_sequence.sv");
}

TEST_P(SVUnitTests, Unsupported_EventType)
{
  expect_encode_throws("event_type.sv");
}

TEST_P(SVUnitTests, Unsupported_RealType)
{
  expect_encode_throws("real_type.sv");
}

TEST_P(SVUnitTests, Unsupported_StringType)
{
  expect_encode_throws("string_type.sv");
}

TEST_P(SVUnitTests, Unsupported_ChandleType)
{
  expect_encode_throws("chandle_type.sv");
}

TEST_P(SVUnitTests, Unsupported_DynamicArray)
{
  expect_encode_throws("dynamic_array.sv");
}

TEST_P(SVUnitTests, Unsupported_QueueType)
{
  expect_encode_throws("queue_type.sv");
}

TEST_P(SVUnitTests, Unsupported_AssocArray)
{
  expect_encode_throws("assoc_array.sv");
}

// This one is a real synthesizable-RTL gap (register files / small
// memories are mainstream, not verification-only), kept Gap_-named even
// though it lives in this file's throw-ledger for organizational
// consistency with the other "throws cleanly" cases.
TEST_P(SVUnitTests, Gap_UnpackedRegfileMemory)
{
  expect_encode_throws("unpacked_regfile.sv");
}

// A register whose output port is aliased through an instance-array
// bus-element connection to only *part* of its target's declared
// width (compare gapped_bus_slice.sv's analogous wire-splicing case)
// isn't supported: declare_variables_internal() has no splicing logic
// for a register spread across sibling instances the way
// process_continuous_assign() does for a wire, so it must throw
// rather than silently never create a state var for the shared target
// at all.
TEST_P(SVUnitTests, Gap_RegisterAliasedToPartialTarget)
{
  expect_encode_throws("reg_bus_slice.sv");
}

// A range-select lvalue with a non-constant (variable) base
// (`w[base +: 4]`) has no dynamic-range-select write fallback anywhere
// in this encoder, unlike ElementSelect's single-bit dynamic-index
// fallback (process_dynamic_element_assign()). resolve_lvalue() throws
// a clear PonoException for this rather than silently dropping the
// write -- the same "throw rather than silently mis-encode" contract
// enforced everywhere else in this file.
TEST_P(SVUnitTests, Gap_DynamicRangeSelectLhs)
{
  expect_encode_throws("dynamic_range_select_lhs.sv");
}

// A streaming concatenation used as an assignment target
// (`{>>{hi, lo}} <= a;`) is ExpressionKind::Streaming, distinct from a
// plain concatenation-target LHS (ExpressionKind::Concatenation,
// already supported). resolve_lvalue() has no case for it at all.
TEST_P(SVUnitTests, Gap_StreamingConcatLhs)
{
  expect_encode_throws("streaming_concat_lhs.sv");
}

// A constant element-select lvalue whose index is out of range for its
// base (`flag[10]` into a 4-bit `flag`) -- the LRM permits this
// (writes are a no-op, reads return 'x), but this encoder has no such
// semantics.
TEST_P(SVUnitTests, Gap_ElementSelectOutOfBoundsLhs)
{
  expect_encode_throws("element_select_out_of_bounds_lhs.sv");
}

// `defparam`/`bind` are legacy simulation-era constructs with no
// functional-logic representation: the base module they target is
// still walked normally with its *own* defaults, so encoding succeeds
// as if the override/bind-in never appeared in the source -- logged
// via logger.log(1, "... ignoring ...") rather than thrown. `defparam`
// is caught as a walkable SymbolKind::DefParam member; `bind` has no
// elaborated Symbol at all and is instead found by scanning the raw
// syntax tree for BindDirective nodes (warn_on_bind_directives()).
TEST_P(SVUnitTests, Unsupported_DefparamStmt)
{
  expect_encode_succeeds_ignoring("defparam_stmt.sv");
}

TEST_P(SVUnitTests, Unsupported_BindDirective)
{
  expect_encode_succeeds_ignoring("bind_directive.sv");
}

// `specify` affects only timing (not functional logic), so ignoring it
// (a walkable SymbolKind::SpecifyBlock member, logged via
// logger.log(1, "... ignoring specify block ...") rather than thrown)
// doesn't corrupt any functional proof the way the assume/cover/
// statement-kind gaps elsewhere in this suite can.
TEST_P(SVUnitTests, Unsupported_SpecifyBlock)
{
  expect_encode_succeeds_ignoring("specify_block.sv");
}

// force/release's effect is inherently simulation-timing-dependent and
// doesn't map onto Pono's per-cycle model; process_statement()'s
// default case logs a warning and skips them (same mechanism as
// fork/join and wait above), so encoding just proceeds as if the
// initial block's force/release calls weren't there (this test's own
// assert doesn't reference `x`, so it can't further distinguish
// "ignored" from "applied and then reverted" -- only that neither one
// crashes the encoder).
TEST_P(SVUnitTests, Unsupported_ForceRelease)
{
  expect_encode_succeeds_ignoring("force_release.sv");
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVUnsupportedTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests

#endif
