#include "sv_test_fixture.h"

using namespace pono;
using namespace smt;

namespace pono_tests {

// ---------------------------------------------------------------------------
// Ledger of genuinely out-of-scope IEEE 1800-2017 constructs: OOP/classes,
// randomization, DPI, functional coverage, programs/specify,
// fork/join/wait/force-release, non-integral types, dynamic containers, and
// hierarchical/port-connection idioms that aren't real synthesizable RTL.
// Every test here is checked via expect_encode_throws() or
// expect_encode_succeeds_ignoring(), whichever matches how the encoder
// actually rejects it -- see the per-test comment when that isn't obvious
// from the test name. The filename already says "unsupported", so test
// names don't repeat it; real synthesizable-RTL/verification-relevant gaps
// (`Gap_`-prefixed tests) live in the topical file matching their subject
// matter instead -- see project_sv_encoder_gaps.md for the full ranked list.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, ClassDecl) { expect_encode_throws("class_decl.sv"); }

TEST_P(SVUnitTests, RandomizeConstraint)
{
  expect_encode_throws("randomize_constraint.sv");
}

TEST_P(SVUnitTests, DpiImport) { expect_encode_throws("dpi_import.sv"); }

TEST_P(SVUnitTests, CovergroupDecl)
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
TEST_P(SVUnitTests, ProgramBlock)
{
  expect_encode_succeeds_ignoring("program_block.sv");
}

// `fork`/`join` and `wait` are simulation-timing constructs with no
// per-cycle counterpart in this encoder's model; process_statement()'s
// default case logs a warning (logger.log(1, "... skipping unsupported
// statement kind ...")) and skips them, rather than throwing.
TEST_P(SVUnitTests, ForkJoin)
{
  expect_encode_succeeds_ignoring("fork_join.sv");
}

TEST_P(SVUnitTests, WaitStmt)
{
  expect_encode_succeeds_ignoring("wait_stmt.sv");
}

// `expect (property_expr);` is a procedural blocking-wait statement
// (pause until the property holds), not a checked invariant -- the
// same simulation-only category as `wait` above.
TEST_P(SVUnitTests, ExpectProperty)
{
  expect_encode_succeeds_ignoring("expect_property.sv");
}

TEST_P(SVUnitTests, EventType) { expect_encode_throws("event_type.sv"); }

TEST_P(SVUnitTests, RealType) { expect_encode_throws("real_type.sv"); }

TEST_P(SVUnitTests, StringType) { expect_encode_throws("string_type.sv"); }

TEST_P(SVUnitTests, ChandleType) { expect_encode_throws("chandle_type.sv"); }

TEST_P(SVUnitTests, DynamicArray) { expect_encode_throws("dynamic_array.sv"); }

TEST_P(SVUnitTests, QueueType) { expect_encode_throws("queue_type.sv"); }

TEST_P(SVUnitTests, AssocArray) { expect_encode_throws("assoc_array.sv"); }

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
TEST_P(SVUnitTests, HierarchicalContinuousAssignForwardRef)
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
TEST_P(SVUnitTests, DynamicIndexOutputPortConnection)
{
  expect_encode_throws("dynamic_index_output_port.sv");
}

// The concatenation-target output-port-connection path
// (`.out({a, bus[idx]})`) must throw as soon as any operand's
// resolve_lvalue() fails, for the same reason as the plain-expression
// case above, rather than silently dropping the whole multi-piece
// write.
TEST_P(SVUnitTests, DynamicIndexConcatOutputPortConnection)
{
  expect_encode_throws("dynamic_index_concat_output_port.sv");
}

// `specify` affects only timing (not functional logic), so ignoring it
// (a walkable SymbolKind::SpecifyBlock member, logged via
// logger.log(1, "... ignoring specify block ...") rather than thrown)
// doesn't corrupt any functional proof the way the assume/cover/
// statement-kind gaps elsewhere in this suite can.
TEST_P(SVUnitTests, SpecifyBlock)
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
TEST_P(SVUnitTests, ForceRelease)
{
  expect_encode_succeeds_ignoring("force_release.sv");
}

// `final` blocks run once at the end of simulation for cleanup/
// reporting -- no synthesis meaning and no analog in this bounded/
// infinite-trace model (there's no "end of simulation"), so they're
// intentionally ignored, the same as $display and other simulation-
// only constructs elsewhere in this encoder.
TEST_P(SVUnitTests, FinalBlockIgnored)
{
  expect_encode_succeeds_ignoring("final_block.sv");
}

// `$display` (and the rest of the display/severity/file-I/O/simulation-
// control system-task families) used as a bare statement inside an
// ordinary procedural block: process_statement()'s Call-expression
// handling recognizes any system call as simulation-only and skips it
// (logged), the same as $display inside the `final` block above --
// distinct from a user-defined task call, which has real side effects
// this encoder can't inline and so throws instead (Gap_VoidTaskCall in
// test_systemverilog_statements.cpp).
TEST_P(SVUnitTests, DisplayCallStatementIgnored)
{
  expect_encode_succeeds_ignoring("display_call_statement.sv");
}

// A bare `forever` (no event control) has no static iteration bound at
// all and can't be unrolled by the compile-time-bounded model -- a
// genuine architectural boundary, not a "not implemented yet" gap. See
// ForeverEventAsRegister in test_systemverilog_statements.cpp for the
// supported `initial forever @(...) ...` structural spelling of a
// register.
TEST_P(SVUnitTests, BareForever) { expect_encode_throws("bare_forever.sv"); }

// `+incdir+`/`-y`-style tool directives in a `.f` list file are
// rejected outright rather than silently treated as filenames --
// parse_dot_f_file() only understands bare filenames and comments,
// a deliberate scope limitation of this encoder's `.f`-file parser
// (an unofficial, multi-vendor EDA build-flow convention, not an IEEE
// 1800 SystemVerilog construct), not a partially-implemented feature.
TEST_P(SVUnitTests, FilelistIncdirDirective)
{
  expect_encode_throws("filelist_top.sv",
                       { sv_path("filelist_bad_directive.f") });
}

// A checker's own formal ports are resolved by slang's elaboration
// (a reference inside the checker body binds directly to the actual
// argument's own symbol, no port-binding work needed on this side --
// see CheckerBlock in test_systemverilog_hierarchy.cpp), but a
// variable declared directly in the checker's own body is genuine new
// local state, needing its own pre-scan/declare pass the way a
// module's does. This encoder doesn't extend those passes into
// checker bodies, so process_checker_instance() throws rather than
// leaving that state undeclared.
TEST_P(SVUnitTests, CheckerLocalVariable)
{
  expect_encode_throws("checker_local_variable.sv");
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVUnsupportedTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
