#ifdef WITH_SLANG

#include "sv_test_fixture.h"

using namespace pono;
using namespace smt;

namespace pono_tests {

TEST_P(SVUnitTests, HierarchicalModules)
{
  check_bmc("hierarchical.sv", 5, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, HierarchicalValue)
{
  check_bmc("hierarchical_value.sv", 6);
}

// An unconnected child input port (`.a()`) has no driver and no
// 4-state representation in this encoder, so it's modeled as a free
// input; `bout`, wired combinationally from it, is freely choosable
// by BMC rather than stuck at a fixed value.
TEST_P(SVUnitTests, UnconnectedInputPortIsFree)
{
  check_bmc("unconnected_input_port.sv", 1);
}

// A child output port connected to a concatenation of parent-side
// signals (`.sum({hi, lo})`) splits the port's bits across two parent
// nets. port_output_aliases_/resolve_output_alias_pieces() map the
// write onto each concatenation operand's own target range, so
// `hi`/`lo` each track their own slice of the child's driven value.
TEST_P(SVUnitTests, ConcatenationOutputPort)
{
  check_bmc("concat_output_port.sv", 4, ProverResult::UNKNOWN);
}

// Three sibling instances drive non-adjacent slices of one shared bus;
// the first-processed instance's write starts at bit 0 without
// covering the whole bus, exercising process_continuous_assign()'s
// full-width-write check against the wire's declared width rather
// than mistaking a partial first write for a full one.
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

TEST_P(SVUnitTests, GenerateIfSelectsFastPath)
{
  check_bmc("generate_if.sv", 4);
}

TEST_P(SVUnitTests, GenerateCaseSelectsBranch)
{
  check_bmc("generate_case.sv", 4);
}

// Parameter override via positional (`#(3)`) vs. named (`#(.WIDTH(3))`)
// syntax.
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

// A plain signal-bundle `interface` wiring a producer module to the
// top level: `bus.data`/`bus.valid` resolve to the interface
// instance's member symbols via the HierarchicalValue path.
// pre_scan_state_vars() recurses into every instance up front, before
// any variables are declared, so the interface's registers (driven by
// a later-declared sibling module) are known regardless of
// declaration order.
TEST_P(SVUnitTests, InterfaceSignalBundle)
{
  check_bmc("interface_bundle.sv", 4);
}

// Same shape as InterfaceSignalBundle, but through a modport-qualified
// port (`bus_if.master b`): `b.data` resolves to a synthesized
// ModportPortSymbol proxy rather than the interface's own `data`
// VariableSymbol directly; canonicalize_modport_port() redirects
// through ModportPortSymbol::internalSymbol so both access paths
// converge on the same underlying symbol.
TEST_P(SVUnitTests, InterfaceModportPort)
{
  check_bmc("interface_modport_task.sv", 4);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVHierarchyTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests

#endif
