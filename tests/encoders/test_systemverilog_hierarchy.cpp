#include "engines/kinduction.h"
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

// A register whose output port is aliased through an instance-array
// bus-element connection to only *part* of its target's declared
// width (compare GappedBusSliceFromSiblingInstances's analogous
// wire-splicing case above) isn't supported:
// declare_variables_internal() has no splicing logic for a register
// spread across sibling instances the way process_continuous_assign()
// does for a wire.
TEST_P(SVUnitTests, Gap_RegisterAliasedToPartialTarget)
{
  check_bmc("reg_bus_slice.sv", 2, ProverResult::UNKNOWN);
}

// A streaming concatenation as an output-port connection. `>>`
// re-orders nothing, so it splits into output-alias segments exactly
// as a plain concatenation connection does.
TEST_P(SVUnitTests, StreamingConcatPortConnection)
{
  check_prover<KInduction>("streaming_concat_port.sv", 4, ProverResult::TRUE);
}

// `<<` moves bits across the boundaries between the stream's
// expressions, and an alias segment is one contiguous range per
// symbol, so there is nothing to describe the result with.
TEST_P(SVUnitTests, Unsupported_StreamingConcatPortReversed)
{
  expect_encode_throws("streaming_concat_port_reversed.sv");
}

TEST_P(SVUnitTests, GenerateForBlock) { check_bmc("generate_block.sv", 6); }

TEST_P(SVUnitTests, ForLoopPopcount) { check_bmc("for_loop.sv", 2); }

TEST_P(SVUnitTests, Parameter) { check_bmc("parameter.sv", 16); }

TEST_P(SVUnitTests, Filelist)
{
  check_bmc(
      "filelist_top.sv", 5, ProverResult::UNKNOWN, { sv_path("filelist.f") });
}

TEST_P(SVUnitTests, FilelistMissingFile)
{
  SmtSolver s = create_solver(GetParam());
  FunctionalTransitionSystem fts(s);
  EXPECT_THROW(
      SystemVerilogEncoder::encode(
          fts, sv_path("filelist_top.sv"), { sv_path("filelist_missing.f") }),
      PonoException);
}

// A module nothing instantiates is top-level, so an unused helper is
// enough to make "the top" ambiguous.  Picking one silently would drop
// the other module -- assertions included -- from the model entirely,
// reporting a clean pass for a design that fails; and it would pick
// badly, since slang orders tops alphabetically rather than by source
// order.  multi_top.sv is named so the alphabetically-first module is
// *not* the intended one.
TEST_P(SVUnitTests, MultipleTopsRejected)
{
  expect_encode_throws("multi_top.sv");
}

// ...and --sv-top resolves it, reaching the module that was previously
// unreachable.  `dut`'s assertion is the false one, so finding exactly
// one property here, and a counterexample for it, is what shows the
// right module got encoded.
TEST_P(SVUnitTests, MultipleTopsSelectedByName)
{
  SmtSolver s = create_solver(GetParam());
  s->set_opt("incremental", "true");
  s->set_opt("produce-models", "true");
  FunctionalTransitionSystem fts(s);
  auto sv_result = SystemVerilogEncoder::encode(
      fts, sv_path("multi_top.sv"), /*filelists=*/{}, /*top=*/"dut");
  ASSERT_EQ(sv_result.propvec.size(), 1u);

  TransitionSystem ts = fts;
  Term prop_term = sv_result.propvec[0];
  Term rst = find_reset(ts);
  ASSERT_TRUE(rst);
  Term reset_done = add_reset_seq(ts, rst, /*reset_bnd=*/1);
  prop_term = ts.solver()->make_term(Implies, reset_done, prop_term);

  SafetyProperty prop(ts.solver(), prop_term);
  Bmc bmc(prop, ts, s);
  // The counter reaches 3 four cycles after reset releases.
  EXPECT_EQ(bmc.check_until(4), ProverResult::FALSE);
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

// `defparam` has real functional effect (it overrides a parameter,
// here changing a counter's bit width) -- not a deliberate non-goal
// like the constructs in test_systemverilog_unsupported.cpp. The base
// module is still walked normally with its *own* defaults, so the
// override is silently never applied -- logged via
// logger.log(1, "... ignoring ...") rather than thrown. `defparam` is
// caught as a walkable SymbolKind::DefParam member.
TEST_P(SVUnitTests, Gap_DefparamStmt)
{
  check_bmc("defparam_stmt.sv", 20, ProverResult::UNKNOWN);
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

// `bind` works: slang's elaborator splices the bound instance directly
// into the target module's body as an ordinary child instance, so
// Pono's regular instance walk picks up `checker_mod`'s assertion --
// and its port connection to leaf2's internal `count` -- for free,
// with no special-case bind handling needed in this encoder at all.
// `warn_on_bind_directives()` (encoder.cpp) only logs an informational
// warning; it doesn't skip anything.
TEST_P(SVUnitTests, BindDirectiveAttachesAssertion)
{
  check_bmc("bind_directive.sv", 3, ProverResult::UNKNOWN);
}

// `checker` is IEEE 1800's standard non-invasive formal-assertion-
// attachment mechanism, not a deliberate non-goal like `program`/
// `specify` in test_systemverilog_unsupported.cpp.
// InstanceEncoder::process_checker_instance() walks a checker
// instance's body exactly like a module's: a reference to one of the
// checker's own formal (AssertionPortSymbol) ports already binds
// directly to the caller's actual argument symbol (slang's own
// elaboration substitutes it), so `a` here resolves straight through
// to checker_block's own free input with no port-binding work needed
// on this side. checker_block has an (unused) `rst` input, so
// check_bmc's own reset-gating pushes the earliest checked cycle to
// 1, not 0.
TEST_P(SVUnitTests, CheckerBlock)
{
  check_bmc("checker_block.sv", 1, ProverResult::FALSE);
}

// A checker instantiated inside a *non-top* module (as opposed to
// checker_block.sv's top-level instantiation) -- exercises
// process_instance()'s own CheckerInstance dispatch, distinct from
// process_assignments()'s top-level one. No `rst` input here, so
// (unlike CheckerBlock above) the earliest violation is at cycle 0.
TEST_P(SVUnitTests, NestedCheckerBlock)
{
  check_bmc("nested_checker_block.sv", 0, ProverResult::FALSE);
}

// A checker with its own local sequential state (`count`, driven by
// an always_ff inside the checker body), not just formal-port
// references: SymbolTable::pre_scan_state_vars() now recurses into
// checker instances to classify `count` as a state var, and
// process_checker_instance() declares it under the checker instance's
// own hierarchical name (checker_local_state.chk.count) exactly like
// a module instance's internal register. After the reset cycle, `a`
// free each cycle can increment `count` from 0 up to 3, violating
// `count < 3`.
TEST_P(SVUnitTests, CheckerLocalState)
{
  check_bmc("checker_local_state.sv", 4, ProverResult::FALSE);
}

// A runtime-indexed write inside a child whose output port is
// connected to a concatenation: the port's bits are spread across
// two parent-side signals, and which one the write reaches is only
// known at runtime, so every segment takes a guarded splice.
TEST_P(SVUnitTests, DynamicWritePortAlias)
{
  check_prover<KInduction>(
      "dynamic_write_port_alias.sv", 8, ProverResult::TRUE);
}

// Asserting only that the right half changed would also pass if
// neither did.
TEST_P(SVUnitTests, DynamicWritePortAliasFails)
{
  check_bmc("dynamic_write_port_alias_fails.sv", 2);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVHierarchyTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

// ---------------------------------------------------------------------------
// `program` instances. A program is a testbench entry point, so its
// stimulus is left unencoded: it drives the DUT along one particular
// scenario, and pinning the inputs to it would leave every other
// input sequence unexplored while still reporting a proof. An
// assertion written inside one is not stimulus, though, and is
// encoded like any other property -- dropping it meant reporting a
// proof over fewer properties than were written.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, ProgramAssertion) { check_bmc("program_assertion.sv", 6); }

TEST_P(SVUnitTests, ProgramAssertionHolds)
{
  check_prover<KInduction>(
      "program_assertion_holds.sv", 12, ProverResult::TRUE);
}

// The stimulus half: this program's body is a `$display` in an
// `initial`, which has no per-cycle meaning and is ignored rather
// than encoded or rejected.
TEST_P(SVUnitTests, ProgramStimulusIgnored)
{
  expect_encode_succeeds_ignoring("program_block.sv");
}

// An `inout` port. Only the child's drive is modelled -- right while
// nothing drives the net from outside, and unable to represent the
// parent driving back, which is logged rather than left silent.

TEST_P(SVUnitTests, InoutPort)
{
  check_prover<KInduction>("inout_port.sv", 8, ProverResult::TRUE);
}

TEST_P(SVUnitTests, InoutPortUndrivenIsFree)
{
  check_bmc("inout_port_fails.sv", 1);
}

}  // namespace pono_tests
