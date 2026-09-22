#include "sv_test_fixture.h"

using namespace pono;
using namespace smt;

namespace pono_tests {

// ---------------------------------------------------------------------------
// Constructs this encoder does not aim to support: non-synthesizable
// SystemVerilog, simulation-only behaviour, and hierarchical idioms
// outside the module-port interface. The filename says "unsupported",
// so test names do not repeat it.
//
// Each asserts how the construct is turned away -- expect_encode_throws()
// where the encoder rejects it, expect_encode_succeeds_ignoring() where
// it is walked past without affecting the model. Which of the two
// applies is a property of the construct rather than an accident:
// something with no hardware meaning at all can be skipped, while
// something that would silently change the model if dropped has to
// throw.
//
// Grouped by why the construct is out of scope; add a new test to the
// group whose reason it shares.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Object-oriented and verification-only language features. None of these
// describe hardware.
// ---------------------------------------------------------------------------

// A class declaration.
TEST_P(SVUnitTests, ClassDecl) { expect_encode_throws("class_decl.sv"); }

// A randomization constraint block.
TEST_P(SVUnitTests, RandomizeConstraint)
{
  expect_encode_throws("randomize_constraint.sv");
}

// A DPI import of a foreign function.
TEST_P(SVUnitTests, DpiImport) { expect_encode_throws("dpi_import.sv"); }

// A covergroup declaration.
TEST_P(SVUnitTests, CovergroupDecl)
{
  expect_encode_throws("covergroup_decl.sv");
}

// ---------------------------------------------------------------------------
// Simulation-time control. These turn on time passing within a clock
// cycle, or on simulation start-up and shutdown, neither of which a
// cycle-based model has. They are walked past rather than rejected:
// with no hardware meaning at all, skipping one cannot change what the
// design does.
// ---------------------------------------------------------------------------

// `fork`/`join` concurrent processes.
TEST_P(SVUnitTests, ForkJoin)
{
  expect_encode_succeeds_ignoring("fork_join.sv");
}

// A bare `wait` statement.
TEST_P(SVUnitTests, WaitStmt)
{
  expect_encode_succeeds_ignoring("wait_stmt.sv");
}

// `expect (property_expr)`, a procedural blocking wait for a property
// rather than a checked invariant.
TEST_P(SVUnitTests, ExpectProperty)
{
  expect_encode_succeeds_ignoring("expect_property.sv");
}

// Procedural `force`/`release`.
TEST_P(SVUnitTests, ForceRelease)
{
  expect_encode_succeeds_ignoring("force_release.sv");
}

// A `final` block, which runs once when simulation ends -- a moment an
// infinite-trace model does not have.
TEST_P(SVUnitTests, FinalBlockIgnored)
{
  expect_encode_succeeds_ignoring("final_block.sv");
}

// `$display` and the rest of the display, severity, file-I/O and
// simulation-control system tasks, called as a bare statement. A call
// to a *user-defined* task is a different matter and throws, since a
// task exists to write its caller's variables.
TEST_P(SVUnitTests, DisplayCallStatementIgnored)
{
  expect_encode_succeeds_ignoring("display_call_statement.sv");
}

// A `specify` block, which carries path delays and timing checks and so
// touches no functional logic.
TEST_P(SVUnitTests, SpecifyBlock)
{
  expect_encode_succeeds_ignoring("specify_block.sv");
}

// ---------------------------------------------------------------------------
// Types with no bit-vector representation: either not integral at all,
// or sized at run time rather than at elaboration.
// ---------------------------------------------------------------------------

// The `event` type.
TEST_P(SVUnitTests, EventType) { expect_encode_throws("event_type.sv"); }

// The `real` type.
TEST_P(SVUnitTests, RealType) { expect_encode_throws("real_type.sv"); }

// The `string` type.
TEST_P(SVUnitTests, StringType) { expect_encode_throws("string_type.sv"); }

// The `chandle` type.
TEST_P(SVUnitTests, ChandleType) { expect_encode_throws("chandle_type.sv"); }

// A dynamic array, whose size is set at run time.
TEST_P(SVUnitTests, DynamicArray) { expect_encode_throws("dynamic_array.sv"); }

// A queue.
TEST_P(SVUnitTests, QueueType) { expect_encode_throws("queue_type.sv"); }

// An associative array.
TEST_P(SVUnitTests, AssocArray) { expect_encode_throws("assoc_array.sv"); }

// A packed struct with a field that is itself an unpacked array: there
// is no flat layout to give that field a bit range in, so the struct
// has no sort. (An unpacked *union* is refused for a different reason
// -- its fields overlap rather than concatenate.)
TEST_P(SVUnitTests, UnpackedStructNonIntegralField)
{
  expect_encode_throws("unpacked_struct_nonintegral.sv");
}

// ---------------------------------------------------------------------------
// Control flow with no static bound. Neither describes a finite piece
// of hardware, and synthesis tools reject them too.
// ---------------------------------------------------------------------------

// A bare `forever` with no event control.
TEST_P(SVUnitTests, BareForever) { expect_encode_throws("bare_forever.sv"); }

// A function that recurses to a depth known only at run time.
// (Recursion to a *constant* depth is a separate question: that
// elaborates statically and could be unrolled.)
TEST_P(SVUnitTests, UnboundedRecursion)
{
  expect_encode_throws("unbounded_recursion.sv");
}

// The same for a task, which is inlined by the same machinery.
TEST_P(SVUnitTests, RecursiveTaskCall)
{
  expect_encode_throws("task_call_recursive.sv");
}

// ---------------------------------------------------------------------------
// Cross-instance wiring outside the module-port interface, which is the
// only sanctioned mechanism; synthesis tools reject these too. A port
// connection is a structural, elaboration-time binding rather than a
// per-cycle write, so there is no mux for a runtime index to build. All
// three throw rather than drop the write, which would leave the target
// unconstrained.
// ---------------------------------------------------------------------------

// A continuous assign driving a child instance's internal, non-port
// signal through a hierarchical dot-path.
TEST_P(SVUnitTests, HierarchicalContinuousAssignForwardRef)
{
  expect_encode_throws("hier_continuous_assign_forward_ref.sv");
}

// An output port connected to a runtime-variable bit select,
// `.out(bus[idx])`.
TEST_P(SVUnitTests, DynamicIndexOutputPortConnection)
{
  expect_encode_throws("dynamic_index_output_port.sv");
}

// The same inside a concatenation target, `.out({a, bus[idx]})`, where
// one unresolvable operand rejects the whole connection.
TEST_P(SVUnitTests, DynamicIndexConcatOutputPortConnection)
{
  expect_encode_throws("dynamic_index_concat_output_port.sv");
}

// ---------------------------------------------------------------------------
// `.f` list files, an unofficial multi-vendor build-flow convention
// rather than an IEEE 1800 construct. The supported syntax is bare
// filenames and comments; every directive is refused rather than
// mistaken for a filename.
// ---------------------------------------------------------------------------

// A `+incdir+`/`-y`-style search-path directive.
TEST_P(SVUnitTests, FilelistIncdirDirective)
{
  expect_encode_throws("filelist_top.sv",
                       { sv_path("filelist_bad_directive.f") });
}

// A `-f other.f` directive nesting one list file inside another.
TEST_P(SVUnitTests, FilelistNestedInclude)
{
  expect_encode_throws("filelist_top.sv", { sv_path("nested_filelist.f") });
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVUnsupportedTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
