#ifdef WITH_SLANG

#include "sv_test_fixture.h"

using namespace pono;
using namespace smt;

namespace pono_tests {

// ---------------------------------------------------------------------------
// Ledger of genuinely out-of-scope IEEE 1800-2017 constructs (OOP/classes,
// randomization, DPI, functional coverage, programs/specify,
// fork/join/wait/force-release, non-integral types, dynamic containers).
// Each is checked via either expect_encode_throws() or
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

// `program` instances are a verification-only construct with no
// functional-logic counterpart: process_instance() recognizes a
// program instance via DefinitionKind::Program and skips it, logged
// via logger.log(1, "... ignoring ... instance ...") rather than
// thrown, per the "simulation-only constructs are dropped and logged"
// half of encode()'s documented contract (see
// SystemVerilogEncoder::encode()'s doc comment).
TEST_P(SVUnitTests, Unsupported_ProgramBlock)
{
  expect_encode_succeeds_ignoring("program_block.sv");
}

// `checker` is IEEE 1800's standard non-invasive formal-assertion-
// attachment mechanism, not a deliberate non-goal like `program` above
// -- a checker instance is a distinct SymbolKind::CheckerInstance the
// usual member walk doesn't match at all, so its own `assert property`
// never reaches the model.
TEST_P(SVUnitTests, Gap_CheckerBlock)
{
  check_bmc("checker_block.sv", 1, ProverResult::FALSE);
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
// it is ever pre-scanned *or* applied, and the skip is logged. A real
// mainstream-RTL gap, not a deliberate non-goal.
TEST_P(SVUnitTests, Gap_PatternCase)
{
  check_bmc("pattern_case.sv", 2, ProverResult::UNKNOWN);
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
// `cover property (a ##1 b)`, which IS supported -- see CoverProperty
// in test_systemverilog_sva.cpp): extending reachability duality
// through the LTL tableau for cover goals is a real gap, not a
// deliberate non-goal. check_bmc() attempts the same reachability
// check CoverProperty uses.
TEST_P(SVUnitTests, Gap_CoverSequence) { check_bmc("cover_sequence.sv", 1); }

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

// A real synthesizable-RTL gap (register files / small memories are
// mainstream, not verification-only): unpacked arrays never build an
// SMT array sort, so the fixture's own read-after-write invariant
// can't be checked.
TEST_P(SVUnitTests, Gap_UnpackedRegfileMemory)
{
  check_bmc("unpacked_regfile.sv", 3, ProverResult::UNKNOWN);
}

// A register whose output port is aliased through an instance-array
// bus-element connection to only *part* of its target's declared
// width (compare gapped_bus_slice.sv's analogous wire-splicing case)
// isn't supported: declare_variables_internal() has no splicing logic
// for a register spread across sibling instances the way
// process_continuous_assign() does for a wire.
TEST_P(SVUnitTests, Gap_RegisterAliasedToPartialTarget)
{
  check_bmc("reg_bus_slice.sv", 2, ProverResult::UNKNOWN);
}

// A range-select lvalue with a non-constant (variable) base
// (`w[base +: 4]`) has no dynamic-range-select write fallback anywhere
// in this encoder, unlike ElementSelect's single-bit dynamic-index
// fallback (process_dynamic_element_assign()). resolve_lvalue() throws
// a clear PonoException for this rather than silently dropping the
// write.
TEST_P(SVUnitTests, Gap_DynamicRangeSelectLhs)
{
  check_bmc("dynamic_range_select_lhs.sv", 2, ProverResult::UNKNOWN);
}

// A streaming concatenation used as an assignment target
// (`{>>{hi, lo}} <= a;`) is ExpressionKind::Streaming, distinct from a
// plain concatenation-target LHS (ExpressionKind::Concatenation,
// already supported). resolve_lvalue() has no case for it at all.
TEST_P(SVUnitTests, Gap_StreamingConcatLhs)
{
  check_bmc("streaming_concat_lhs.sv", 2, ProverResult::UNKNOWN);
}

// A constant element-select lvalue whose index is out of range for its
// base (`flag[10]` into a 4-bit `flag`) -- the LRM permits this
// (writes are a no-op, reads return 'x), but this encoder has no such
// semantics.
TEST_P(SVUnitTests, Gap_ElementSelectOutOfBoundsLhs)
{
  check_bmc("element_select_out_of_bounds_lhs.sv", 2, ProverResult::UNKNOWN);
}

// A continuous assign targeting a child instance's internal (non-port)
// signal via a hierarchical dot-path is not real synthesizable RTL to
// begin with -- module ports are the only sanctioned cross-instance
// wiring mechanism, so driving a submodule's internals directly from
// outside its own scope is a simulation/testbench/hierarchical-deposit
// idiom, a deliberate non-goal rather than missed synthesizable-subset
// work. This particular fixture also appears before that instance's
// own declaration in the same scope (declaration is interleaved with,
// and ordered by, source position), so the target has no declared term
// yet when process_continuous_assign_operand() processes it -- it
// throws rather than silently dropping the write either way.
TEST_P(SVUnitTests, Unsupported_HierarchicalContinuousAssignForwardRef)
{
  expect_encode_throws("hier_continuous_assign_forward_ref.sv");
}

// An output port connected via a dynamic (runtime-variable) bit-select
// (`.out(bus[idx])`) isn't real synthesizable RTL to begin with: a
// port connection is a structural, elaboration-time binding, not a
// per-cycle write, so there's no mux for a dynamic index to build --
// a deliberate non-goal rather than missed synthesizable-subset work,
// like the hierarchical-reference case above. resolve_lvalue() throws
// a clear PonoException rather than silently dropping the child's
// output write and leaving the target fully unconstrained.
TEST_P(SVUnitTests, Unsupported_DynamicIndexOutputPortConnection)
{
  expect_encode_throws("dynamic_index_output_port.sv");
}

// The concatenation-target output-port-connection path
// (`.out({a, bus[idx]})`) must throw as soon as any operand's
// resolve_lvalue() fails, for the same reason as the plain-expression
// case above, rather than silently dropping the whole multi-piece
// write.
TEST_P(SVUnitTests, Unsupported_DynamicIndexConcatOutputPortConnection)
{
  expect_encode_throws("dynamic_index_concat_output_port.sv");
}

// `defparam` has real functional effect (it overrides a parameter,
// here changing a counter's bit width) -- not a deliberate non-goal
// like the constructs above. The base module is still walked normally
// with its *own* defaults, so the override is silently never applied
// -- logged via logger.log(1, "... ignoring ...") rather than thrown.
// `defparam` is caught as a walkable SymbolKind::DefParam member.
TEST_P(SVUnitTests, Gap_DefparamStmt)
{
  check_bmc("defparam_stmt.sv", 20, ProverResult::UNKNOWN);
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

// `final` blocks run once at the end of simulation for cleanup/
// reporting -- no synthesis meaning and no analog in this bounded/
// infinite-trace model (there's no "end of simulation"), so they're
// intentionally ignored, the same as $display and other simulation-
// only constructs elsewhere in this encoder.
TEST_P(SVUnitTests, Unsupported_FinalBlockIgnored)
{
  expect_encode_succeeds_ignoring("final_block.sv");
}

// A bare `forever` (no event control) has no static iteration bound at
// all and can't be unrolled by the compile-time-bounded model -- a
// genuine architectural boundary, not a "not implemented yet" gap. See
// ForeverEventAsRegister in test_systemverilog_statements.cpp for the
// supported `initial forever @(...) ...` structural spelling of a
// register.
TEST_P(SVUnitTests, Unsupported_BareForever)
{
  expect_encode_throws("bare_forever.sv");
}

// `+incdir+`/`-y`-style tool directives in a `.f` list file are
// rejected outright rather than silently treated as filenames --
// parse_dot_f_file() only understands bare filenames and comments,
// a deliberate scope limitation of this encoder's `.f`-file parser
// (an unofficial, multi-vendor EDA build-flow convention, not an IEEE
// 1800 SystemVerilog construct), not a partially-implemented feature.
TEST_P(SVUnitTests, Unsupported_FilelistIncdirDirective)
{
  expect_encode_throws("filelist_top.sv",
                       { sv_path("filelist_bad_directive.f") });
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVUnsupportedTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests

#endif
