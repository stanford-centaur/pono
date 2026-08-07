#ifdef WITH_SLANG

#include "sv_test_fixture.h"

using namespace pono;
using namespace smt;

namespace pono_tests {

// ---------------------------------------------------------------------------
// Ledger of genuinely out-of-scope IEEE 1800-2017 chapters (OOP/classes,
// randomization, DPI, functional coverage, programs/checkers,
// fork/join/wait, non-integral types, dynamic containers) plus a few
// synthesizable-RTL legacy constructs (defparam, bind).  Each is checked
// via expect_encode_throws() per the intended "clean rejection" contract;
// where the actual behavior turned out to be something else (a silent
// drop, or a different exception message), that's called out in the
// per-test comment -- verified empirically, not assumed.  All test names
// are prefixed `Unsupported_` (rather than `Gap_`) since these are
// deliberate non-goals, not missed synthesizable-subset features -- with
// the exception of unpacked-array memories and plain interfaces, which
// stayed in Gap_ ledgers (test_systemverilog_unsupported.cpp /
// test_systemverilog_hierarchy.cpp respectively) since those *are*
// mainstream RTL, not verification-only features.
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

// GAP in the "clean rejection" contract itself (confirmed empirically):
// program/checker instances are never walked at all -- process_module()'s
// top-level walk and process_instance() only recognize ordinary module
// InstanceSymbols -- so the whole construct is silently invisible rather
// than throwing.
TEST_P(SVUnitTests, Unsupported_ProgramBlock)
{
  expect_encode_succeeds_ignoring("program_block.sv");
}

TEST_P(SVUnitTests, Unsupported_CheckerBlock)
{
  expect_encode_succeeds_ignoring("checker_block.sv");
}

// GAP in the "clean rejection" contract (confirmed empirically): these
// hit the same silent statement-kind skip documented in
// test_systemverilog_statements.cpp, not a throw.
TEST_P(SVUnitTests, Unsupported_ForkJoin)
{
  expect_encode_succeeds_ignoring("fork_join.sv");
}

TEST_P(SVUnitTests, Unsupported_WaitStmt)
{
  expect_encode_succeeds_ignoring("wait_stmt.sv");
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

// GAP in the "clean rejection" contract (confirmed empirically):
// `defparam`/`bind` are simply never processed -- the base module they
// target is still walked normally with its *own* defaults, so encoding
// succeeds as if the override/bind-in never appeared in the source.
TEST_P(SVUnitTests, Unsupported_DefparamStmt)
{
  expect_encode_succeeds_ignoring("defparam_stmt.sv");
}

TEST_P(SVUnitTests, Unsupported_BindDirective)
{
  expect_encode_succeeds_ignoring("bind_directive.sv");
}

// GAP in the "clean rejection" contract (confirmed empirically, and
// distinct from the plain signal-bundle interface case in
// test_systemverilog_hierarchy.cpp, which *does* throw "Unknown state
// variable" building the same shape of property): with a modport-
// qualified port connection (`b.master`), encoding succeeds without
// throwing at all, even though the assert references `b.data` directly
// -- i.e. this shows a *different* failure mode than the plain-bundle
// case, not just a more severe version of it.
TEST_P(SVUnitTests, Unsupported_InterfaceModportTask)
{
  expect_encode_succeeds_ignoring("interface_modport_task.sv");
}

// `specify` affects only timing (not functional logic), so silently
// ignoring it (confirmed empirically: it's a module-scope member kind
// walk_members() doesn't recognize, skipped with no error) doesn't
// corrupt any functional proof the way the assume/cover/statement-kind
// gaps elsewhere in this suite can.
TEST_P(SVUnitTests, Unsupported_SpecifyBlock)
{
  expect_encode_succeeds_ignoring("specify_block.sv");
}

// force/release's effect is inherently simulation-timing-dependent and
// doesn't map onto Pono's per-cycle model; confirmed empirically that
// encoding just proceeds as if the initial block's force/release calls
// weren't there (this test's own assert doesn't reference `x`, so it
// can't further distinguish "ignored" from "applied and then reverted" --
// only that neither one crashes the encoder).
TEST_P(SVUnitTests, Unsupported_ForceRelease)
{
  expect_encode_succeeds_ignoring("force_release.sv");
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVUnsupportedTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests

#endif
