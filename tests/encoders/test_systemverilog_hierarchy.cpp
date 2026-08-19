#ifdef WITH_SLANG

#include "sv_test_fixture.h"

using namespace pono;
using namespace smt;

namespace pono_tests {

// ---------------------------------------------------------------------------
// Migrated from the original test_systemverilog.cpp -- already compositional,
// kept as-is per the triage in the plan.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, HierarchicalModules)
{
  check_bmc("hierarchical.sv", 5, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, HierarchicalValue)
{
  check_bmc("hierarchical_value.sv", 6);
}

// A child instance's input port left explicitly unconnected (`.a()`):
// declare_variables_internal()'s generic per-body Variable handling
// falls back to a free/undriven input for any symbol that isn't
// classified as a wire, state variable, or output-port alias -- which
// is exactly what an unconnected input port's internal symbol is, so
// it's already treated as a free input rather than crashing or
// silently defaulting to a fixed value. `bout` (wired combinationally
// from the child's now-free `a`) is freely choosable by BMC.
TEST_P(SVUnitTests, UnconnectedInputPortIsFree)
{
  check_bmc("unconnected_input_port.sv", 1);
}

// GAP: a child instance's output port connected to a concatenation of
// parent-side signals (`.sum({hi, lo})`), splitting the port's bits
// across two parent nets. find_lhs_base()/resolve_lvalue() only
// recognize a plain name/index/range/member as a port-connection
// target; a Concatenation expression falls through both (returns
// nullptr/nullopt), so the whole connection is silently dropped --
// `hi`/`lo` end up free rather than tracking the child's driven value,
// despite a comment at the call site already (incorrectly) claiming
// this case "is not yet supported" as if it threw.
TEST_P(SVUnitTests, Gap_ConcatenationOutputPort)
{
  check_bmc("concat_output_port.sv", 4, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, GenerateForBlock) { check_bmc("generate_block.sv", 6); }

TEST_P(SVUnitTests, ForLoopPopcount) { check_bmc("for_loop.sv", 2); }

TEST_P(SVUnitTests, Parameter) { check_bmc("parameter.sv", 16); }

TEST_P(SVUnitTests, Filelist)
{
  check_bmc(
      "filelist_top.sv", 5, ProverResult::UNKNOWN, { sv_path("filelist.f") });
}

TEST_P(SVUnitTests, FilelistUnsupportedDirective)
{
  SmtSolver s = create_solver(GetParam());
  FunctionalTransitionSystem fts(s);
  EXPECT_THROW(
      SystemVerilogEncoder enc(sv_path("filelist_top.sv"),
                               fts,
                               { sv_path("filelist_bad_directive.f") }),
      PonoException);
}

TEST_P(SVUnitTests, FilelistMissingFile)
{
  SmtSolver s = create_solver(GetParam());
  FunctionalTransitionSystem fts(s);
  EXPECT_THROW(
      SystemVerilogEncoder enc(
          sv_path("filelist_top.sv"), fts, { sv_path("filelist_missing.f") }),
      PonoException);
}

// ---------------------------------------------------------------------------
// Generate-if / generate-case: only generate-for was previously tested.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, GenerateIfSelectsFastPath)
{
  check_bmc("generate_if.sv", 4);
}

TEST_P(SVUnitTests, GenerateCaseSelectsBranch)
{
  check_bmc("generate_case.sv", 4);
}

// ---------------------------------------------------------------------------
// Positional vs. named parameter override -- the pre-existing
// parameter.sv test only ever used the module's *default* parameter
// value, never actually overriding it.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, ParameterPositionalOverride)
{
  check_bmc("positional_param_override.sv", 8);
}

TEST_P(SVUnitTests, ParameterNamedOverride)
{
  check_bmc("named_param_override.sv", 8);
}

// ---------------------------------------------------------------------------
// Compiler macros (`define`/`ifdef`) and package `import`, both resolved
// across a --sv-filelist boundary.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, PackageImportMacrosAcrossFilelist)
{
  check_bmc("pkg_macro_top.sv",
            4,
            ProverResult::FALSE,
            { sv_path("pkg_filelist.f") });
}

// ---------------------------------------------------------------------------
// $clog2/$bits as elaboration-time functions sizing an RTL register.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, ClogTwoAndBitsElaborationFunctions)
{
  check_bmc("elaboration_functions.sv", 5);
}

// ---------------------------------------------------------------------------
// A plain signal-bundle `interface` wiring a producer module to the top
// level: `bus.data`/`bus.valid` (referenced from `bus_producer`'s
// always_ff) resolve to the interface instance's real member symbols
// via the HierarchicalValue path. pre_scan_state_vars() recursively
// classifies every always_ff/always block anywhere in the whole
// instance tree *before* any instance's variables are declared, so a
// sibling instance's registers (here, the interface's own members,
// driven by a later-declared sibling module) are known regardless of
// declaration order.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, InterfaceSignalBundle)
{
  check_bmc("interface_bundle.sv", 4);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVHierarchyTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests

#endif
