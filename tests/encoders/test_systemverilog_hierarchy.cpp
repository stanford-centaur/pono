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

// A child instance's output port connected to a concatenation of
// parent-side signals (`.sum({hi, lo})`), splitting the port's bits
// across two parent nets. port_output_aliases_ maps each port-internal
// symbol to one or more OutputAliasSegments (one per concatenation
// operand, MSB-first); resolve_output_alias_pieces() intersects a
// write's own bit window against each segment and remaps it into the
// segment's own target range, so `hi`/`lo` each track their own slice
// of the child's driven value.
TEST_P(SVUnitTests, ConcatenationOutputPort)
{
  check_bmc("concat_output_port.sv", 4, ProverResult::UNKNOWN);
}

// Three sibling instances drive non-adjacent slices of one shared bus,
// with the first-processed instance's write starting at bit 0 but not
// covering the whole bus. process_continuous_assign() must seed a
// full-width placeholder for the wire's first write in that case,
// rather than mistaking it for a full-width write and later corrupting
// (or crashing while building) a subsequent non-adjacent slice write.
TEST_P(SVUnitTests, GappedBusSliceFromSiblingInstances)
{
  check_bmc("gapped_bus_slice.sv", 4, ProverResult::UNKNOWN);
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

// Same shape as InterfaceSignalBundle, but through a *modport*-qualified
// port (`bus_if.master b`) instead of a plain interface port. A
// modport-qualified access (`b.data`) resolves to a synthesized
// ModportPortSymbol proxy rather than directly to the interface
// instance's own `data` VariableSymbol the way the plain (non-modport)
// case does -- canonicalize_modport_port() redirects through
// ModportPortSymbol::internalSymbol so both access paths converge on
// the same underlying symbol identity.
TEST_P(SVUnitTests, InterfaceModportPort)
{
  check_bmc("interface_modport_task.sv", 4);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVHierarchyTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests

#endif
