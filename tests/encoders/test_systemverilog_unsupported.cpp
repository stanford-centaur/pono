#include "engines/kinduction.h"
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
// Each asserts how the construct is turned away: expect_encode_throws()
// where the encoder rejects it, and where it is walked past instead, a
// property that is still proved (and a companion that is still
// refuted) with the construct in scope, which is what shows it made no
// difference. Which of the two applies is a property of the construct
// rather than an accident: something with no hardware meaning at all
// can be skipped, while something that would silently change the model
// if dropped has to throw. A few older walked-past tests still assert
// only that encoding succeeded.
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

// A covergroup, the one construct in this group that is walked past
// rather than rejected: functional coverage observes the design and
// never drives it, so nothing the model keeps depends on it.
TEST_P(SVUnitTests, CovergroupIgnored)
{
  check_prover<KInduction>("covergroup_ignored.sv", 14, ProverResult::TRUE);
}

// The same design and covergroup under a bound the counter passes, so
// the covergroup is shown not to suppress a counterexample either.
TEST_P(SVUnitTests, CovergroupIgnoredBoundFails)
{
  check_bmc("covergroup_ignored_bound.sv", 7);
}

// ---------------------------------------------------------------------------
// Simulation-time control. These turn on time passing within a clock
// cycle, or on simulation start-up and shutdown, neither of which a
// cycle-based model has. Which of them can be walked past divides on
// whether anything the model keeps depends on the timing: a `final`
// block or a `$display` can be dropped whole, while a `fork` or a
// `wait` decides when the statements around it run, and dropping one
// answers that question by accident.
// ---------------------------------------------------------------------------

// `fork`/`join` and its two other terminators. Each statement in a
// parallel block is a concurrent process with timing of its own,
// which a cycle-based transition system has no notion of. Running
// the branches one after another instead would decide a race the
// language leaves open, and for `join_none` would run code that has
// not started.
TEST_P(SVUnitTests, ForkJoin) { expect_encode_throws("fork_join.sv"); }

TEST_P(SVUnitTests, ForkJoinAny) { expect_encode_throws("fork_join_any.sv"); }

TEST_P(SVUnitTests, ForkJoinNone) { expect_encode_throws("fork_join_none.sv"); }

// The shape the construct exists for: branches whose delays are not
// ordered against each other.
TEST_P(SVUnitTests, ForkJoinTimed)
{
  expect_encode_throws("fork_join_timed.sv");
}

// `wait (expr)`, level-sensitive event control: the statements it
// guards stay blocked until the condition holds. Parking a process
// part-way through and resuming it later is per-process state that a
// transition system does not carry. The two spellings differ only in
// whether the guarded write is the wait's own substatement or sits
// after it in the enclosing block; they mean the same thing.
TEST_P(SVUnitTests, WaitLevel) { expect_encode_throws("wait_level.sv"); }

TEST_P(SVUnitTests, WaitThen) { expect_encode_throws("wait_then.sv"); }

// The level-sensitive latch idiom, the one shape of `wait` that
// describes hardware rather than a testbench.
TEST_P(SVUnitTests, WaitLatch) { expect_encode_throws("wait_latch.sv"); }

// `wait fork`, which waits on child processes rather than on a value.
TEST_P(SVUnitTests, WaitFork) { expect_encode_throws("wait_fork.sv"); }

// `wait_order`, which waits for named events to arrive in a given
// order. Its fixture's module cannot be called `wait_order`, since
// that is a keyword.
TEST_P(SVUnitTests, WaitOrder) { expect_encode_throws("wait_order_stmt.sv"); }

// `expect (property_expr)`, a procedural blocking wait for a property
// rather than a checked invariant.
TEST_P(SVUnitTests, ExpectProperty)
{
  expect_encode_succeeds_ignoring("expect_property.sv");
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
// or sized at run time rather than at elaboration. Declaring a
// variable of one is allowed and computes nothing -- see
// UnreferencedUnrepresentableDecls in test_systemverilog_types.cpp --
// so each fixture below goes on to reference its variable, which is
// the first point at which a value is actually needed.
// ---------------------------------------------------------------------------

// Triggering an `event` with `->`.
TEST_P(SVUnitTests, EventType) { expect_encode_throws("event_type.sv"); }

// Reading a `real` variable.
TEST_P(SVUnitTests, RealType) { expect_encode_throws("real_type.sv"); }

// Writing a `string` variable.
TEST_P(SVUnitTests, StringType) { expect_encode_throws("string_type.sv"); }

// Writing a `chandle`.
TEST_P(SVUnitTests, ChandleType) { expect_encode_throws("chandle_type.sv"); }

// A dynamic array, whose size is set at run time.
TEST_P(SVUnitTests, DynamicArray) { expect_encode_throws("dynamic_array.sv"); }

// A queue.
TEST_P(SVUnitTests, QueueType) { expect_encode_throws("queue_type.sv"); }

// An associative array.
TEST_P(SVUnitTests, AssocArray) { expect_encode_throws("assoc_array.sv"); }

// ---------------------------------------------------------------------------
// Control flow with no static bound. Neither describes a finite piece
// of hardware, and synthesis tools reject them too.
// ---------------------------------------------------------------------------

// A bare `forever` with no event control.
TEST_P(SVUnitTests, BareForever) { expect_encode_throws("bare_forever.sv"); }

// A `force` held for a delay rather than for a number of cycles.
// The clock-scoped spellings of force and release are modelled.
TEST_P(SVUnitTests, ForceReleaseWithDelay)
{
  expect_encode_throws("force_release_delay.sv");
}

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
// A runtime index in a port connection. A port connection is a
// structural, elaboration-time binding rather than a per-cycle write,
// so there is no mux for an index to build, and synthesis tools
// reject these too. Both throw rather than drop the write, which
// would leave the target unconstrained.
// ---------------------------------------------------------------------------

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
