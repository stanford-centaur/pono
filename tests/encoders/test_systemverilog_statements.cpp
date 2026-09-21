#include "engines/kinduction.h"
#include "sv_test_fixture.h"

using namespace pono;
using namespace smt;

namespace pono_tests {

// ---------------------------------------------------------------------------
// Compositional cases: a procedural for-loop combined with compound
// assignment (`|=`) or an element-select LHS, both inside always_ff.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, CompoundAssignOrReduce)
{
  check_bmc("compound_assign.sv", 2);
}

TEST_P(SVUnitTests, ElementSelectLhs) { check_bmc("element_select_lhs.sv", 2); }

// Range-select (bit-slice) LHS on a plain continuous assign
// (`assign w[7:4] = ...;`): resolve_lvalue() has a dedicated
// RangeSelect case requiring constant bounds, so each assign pins its
// half of `w` instead of leaving it a free state var.
TEST_P(SVUnitTests, RangeSelectLhs)
{
  check_bmc("range_select_lhs.sv", 4, ProverResult::UNKNOWN);
}

// A range-select lvalue with a non-constant base (`w[base +: 4]`):
// a fixed-width window at a runtime position, which is the same
// splice a runtime-indexed element select needs, only counted in
// bits rather than elements. The read side shifts the window down
// and truncates.
TEST_P(SVUnitTests, DynamicRangeSelectLhs)
{
  check_prover<KInduction>(
      "dynamic_range_select_lhs.sv", 12, ProverResult::TRUE);
}

// A constant element-select lvalue whose index is out of range for its
// base (`flag[10]` into a 4-bit `flag`). The LRM makes the write a
// no-op, which the dynamic-position splice gives for free: a position
// past the end shifts the write mask away entirely.
TEST_P(SVUnitTests, ElementSelectOutOfBoundsLhs)
{
  check_prover<KInduction>(
      "element_select_out_of_bounds_lhs.sv", 6, ProverResult::TRUE);
}

// A runtime-indexed write into a packed range that does not start at
// zero (`logic [7:4] r; r[i] <= ...`). The declared index is a bit
// position only for an `[n:0]` range, so without rebasing the write
// lands four bits too high and quietly misses.
TEST_P(SVUnitTests, DynamicWriteNonZeroBasedRange)
{
  check_prover<KInduction>("dynamic_write_rebased.sv", 8, ProverResult::TRUE);
}

// The ascending twin, where the rebasing runs the other way: for
// `logic [4:7] r`, bit position and declared index count in opposite
// directions.
TEST_P(SVUnitTests, DynamicWriteAscendingRange)
{
  check_prover<KInduction>("dynamic_write_ascending.sv", 8, ProverResult::TRUE);
}

// The read side of the same thing, which the write path's guard never
// covered: bits the vector does not have read as X. Whatever the
// select does reach must be untouched, which is what this proves --
// the two refutations below are what stop "X" being modelled as some
// fixed value.
TEST_P(SVUnitTests, SelectOutOfRangeReads)
{
  check_prover<KInduction>("select_out_of_range.sv", 6, ProverResult::TRUE);
}

TEST_P(SVUnitTests, SelectOutOfRangeConstantFails)
{
  check_bmc("select_out_of_range_fails.sv", 0);
}

// This one was provable before: the shift fed zeros in past the end.
TEST_P(SVUnitTests, SelectOutOfRangeDynamicFails)
{
  check_bmc("select_out_of_range_dynamic_fails.sv", 0);
}

// Concatenation-target LHS on a plain continuous assign (`assign {hi,
// lo} = ...;`), as opposed to a concatenation-target *port connection*
// (already supported separately via OutputAliasSegment). Since a
// concatenation has more than one base symbol and can't be represented
// as a single LValueDesc, process_continuous_assign() special-cases it
// by splitting the RHS across each operand rather than going through
// resolve_lvalue() directly.
TEST_P(SVUnitTests, ConcatenationLhs)
{
  check_bmc("concat_lhs.sv", 4, ProverResult::UNKNOWN);
}

// The wire pre-scan (process_module(), encoder.cpp) that decides
// whether a continuous-assign LHS's base symbol is a wire has its own,
// separate LHS-classification logic from begin_write()'s write-time
// handling above -- confirming `hi`/`lo` are classified as wires
// (macro-substituted, absent from inputvars()) rather than silently
// falling through to free input vars is not otherwise observable via
// check_bmc(), since the write-processing fallback path for an
// unclassified variable still constrains it correctly via
// add_constraint(); only the resulting inputvars()/named_terms()
// bookkeeping differs.
TEST_P(SVUnitTests, ConcatenationLhsClassifiedAsWire)
{
  SmtSolver s = create_solver(GetParam());
  FunctionalTransitionSystem fts(s);
  SystemVerilogEncoder::encode(fts, sv_path("concat_lhs.sv"));
  TransitionSystem ts = fts;
  const auto & named = ts.named_terms();
  ASSERT_TRUE(named.count("concat_lhs.hi"));
  ASSERT_TRUE(named.count("concat_lhs.lo"));
  ASSERT_TRUE(named.count("concat_lhs.a"));
  ASSERT_TRUE(named.count("concat_lhs.b"));
  EXPECT_FALSE(ts.inputvars().count(named.at("concat_lhs.hi")));
  EXPECT_FALSE(ts.inputvars().count(named.at("concat_lhs.lo")));
  EXPECT_TRUE(ts.inputvars().count(named.at("concat_lhs.a")));
  EXPECT_TRUE(ts.inputvars().count(named.at("concat_lhs.b")));
}

// Same shape, procedural (non-blocking) assignment form: begin_write()
// (shared by blocking/non-blocking assignment and ++/--) special-cases
// a top-level concatenation-target LHS the same way
// process_continuous_assign() does for the continuous-assign form
// above.
TEST_P(SVUnitTests, ConcatenationLhsNextState)
{
  check_bmc("concat_lhs_next_state.sv", 2);
}

// A register pair written ONLY through a concat-target NB assignment
// (no reset branch, no other plain-assignment write path anywhere
// else) -- unlike concat_lhs_next_state.sv above, where hi/lo are also
// written via a plain (non-concat) NB assign in the reset branch, so
// they get classified as state vars through that unrelated path
// regardless. This isolates the pre-scan classification gap: the
// concat-target write's own base symbols must be recognized as state
// vars by collect_nonblocking_targets() itself.
TEST_P(SVUnitTests, ConcatenationLhsOnlyWrite)
{
  check_bmc("concat_lhs_only_write.sv", 1);
}

// A streaming concatenation used as an assignment target
// (`{>>{hi, lo}} <= a;`) is ExpressionKind::Streaming, distinct from
// a plain concatenation-target LHS (ExpressionKind::Concatenation,
// handled above) -- but it splits the same way, since `>>` re-orders
// nothing.
TEST_P(SVUnitTests, StreamingConcatLhs)
{
  check_prover<KInduction>("streaming_concat_lhs.sv", 6, ProverResult::TRUE);
}

// The `<<` direction, with a slice size that does not divide the
// width -- where packing and unpacking are different permutations
// rather than the same one applied twice -- and a target streaming
// two expressions rather than one.
TEST_P(SVUnitTests, StreamingConcatReorder)
{
  check_prover<KInduction>(
      "streaming_concat_reorder.sv", 6, ProverResult::TRUE);
}

// The same target in a continuous assignment, with a source wider
// than the target: a stream is consumed from its most significant
// end, not truncated at the bottom.
TEST_P(SVUnitTests, StreamingConcatContinuousAssign)
{
  check_prover<KInduction>(
      "streaming_concat_continuous.sv", 4, ProverResult::TRUE);
}

// Streaming an unpacked array means walking its elements in `foreach`
// order, which positional bit-splicing cannot express.
TEST_P(SVUnitTests, Unsupported_StreamingConcatUnpackedTarget)
{
  expect_encode_throws("streaming_concat_unpacked_target.sv");
}

// A `with` range sizes dynamically sized stream data, which is out of
// scope here.
TEST_P(SVUnitTests, Unsupported_StreamingConcatWithRange)
{
  expect_encode_throws("streaming_concat_with_range.sv");
}

// Minimal, direct checks of two patterns that recur composed with other
// constructs throughout this suite: a bare always_ff counter, and (for
// initial_block.sv below) a design with no `rst` port at all, relying
// solely on `initial` to pin state.
TEST_P(SVUnitTests, EncodeCounter) { check_bmc("counter.sv", 5); }

TEST_P(SVUnitTests, InitialBlockSetsState) { check_bmc("initial_block.sv", 0); }

// `always_latch`: pre_scan_always_latch() marks every blocking-
// assignment target as a state variable unconditionally (unlike
// always_comb's full-vs-partial wire/state-var split), and
// process_next_state_body() (shared with always_ff) processes the body
// with StmtContext::NEXT_STATE, so an unwritten path implicitly holds
// the latch's previous value -- the same "defaults to itself"
// semantics a register's next-state gets.
TEST_P(SVUnitTests, AlwaysLatchHold)
{
  check_bmc("always_latch.sv", 4, ProverResult::UNKNOWN);
}

// A blocking write in a clocked block infers a register too, and is
// visible to later reads in that same block. Proving this rather than
// failing to refute it is what separates the LRM's semantics from the
// non-blocking reading, which would refute every conjunct but the
// self-update.
TEST_P(SVUnitTests, BlockingAssignInClockedBlock)
{
  check_prover<KInduction>(
      "blocking_in_clocked_block.sv", 12, ProverResult::TRUE);
}

// Only its event control says whether a plain `always` writing with
// `=` is a flop or combinational logic; the fixture asserts both
// readings at once, so treating either as the other is refuted.
TEST_P(SVUnitTests, BlockingAlwaysEdgeVsLevel)
{
  check_prover<KInduction>("blocking_always_edge.sv", 12, ProverResult::TRUE);
}

// Writes to one non-wire always_comb target compose in order rather
// than each constraining its own slice. Separate constraints would
// contradict each other and make the design vacuous, which is why
// this is proved and paired with a refutation below.
TEST_P(SVUnitTests, CombPartialWriteComposes)
{
  check_prover<KInduction>("comb_partial_write.sv", 8, ProverResult::TRUE);
}

TEST_P(SVUnitTests, CombPartialWriteFails)
{
  check_bmc("comb_partial_write_fails.sv", 0);
}

// An `initial` block composes its writes the same way, and for the
// same reason: separate per-write constraints leave no satisfiable
// initial state, which refutes nothing at all.
TEST_P(SVUnitTests, InitialPartialWriteComposes)
{
  check_prover<KInduction>("initial_partial_write.sv", 6, ProverResult::TRUE);
}

TEST_P(SVUnitTests, InitialPartialWriteFails)
{
  check_bmc("initial_partial_write_fails.sv", 0);
}

// `p[2][j] <= val`: the base of the runtime-indexed write is itself a
// select, so the splice happens at that inner range's offset. The
// untouched neighbours are asserted too -- an offset error would move
// the write into one of them.
TEST_P(SVUnitTests, NestedDynamicBitWrite)
{
  check_prover<KInduction>(
      "nested_dynamic_bit_write.sv", 12, ProverResult::TRUE);
}

// Both of these used to drop the write with no diagnostic, leaving
// the target free.
TEST_P(SVUnitTests, ConcatDynamicOperandRejected)
{
  expect_encode_throws("concat_dynamic_operand.sv");
}

TEST_P(SVUnitTests, InitialDynamicIndexRejected)
{
  expect_encode_throws("initial_dynamic_index.sv");
}

// A temporary declared inside a procedural block is bound to the
// value it carries, not turned into a register. An uninitialized
// 4-state one used to be forced to its all-X default and stringified
// into the solver, which aborted with no source location at all.
TEST_P(SVUnitTests, ProceduralTemporary)
{
  check_prover<KInduction>("procedural_temp.sv", 12, ProverResult::TRUE);
}

// The other half: a local read where no path assigned it holds its
// value instead, which is storage rather than an unknown. Proving
// this is what separates holding from a fresh unconstrained value
// each cycle, which would also encode but admit traces the hardware
// cannot produce.
TEST_P(SVUnitTests, ProceduralTemporaryHolds)
{
  check_prover<KInduction>("procedural_temp_hold.sv", 12, ProverResult::TRUE);
}

// `initial forever @(posedge clk) ...` is a legacy structural spelling
// of `always_ff @(posedge clk) ...`: as_forever_event_body() recognizes
// this shape (a ForeverLoop whose own body is a Timed statement) and
// redirects it to the same process_next_state_body() an always_ff
// block gets, instead of treating it as an initial-state constraint. A
// *bare* `forever` (no event control) doesn't match this shape and
// remains an architectural boundary -- see BareForever in
// test_systemverilog_unsupported.cpp.
TEST_P(SVUnitTests, ForeverEventAsRegister)
{
  check_bmc("forever_loop.sv", 4, ProverResult::UNKNOWN);
}

// ---------------------------------------------------------------------------
// `priority if` / `unique case` as semantic modifiers (not just a parse-
// through of plain if/case).
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, PriorityIfUniqueCaseFails)
{
  check_bmc("priority_if_unique_case.sv", 2);
}

TEST_P(SVUnitTests, PriorityIfUniqueCaseHolds)
{
  check_bmc("priority_if_unique_case_holds.sv", 6, ProverResult::UNKNOWN);
}

// A `&&&`-joined multi-condition `if` (`if (a &&& b) ...`, LRM 12.4.4)
// must AND together every joined condition rather than reading only
// conditions[0]; mirrors MultiConditionTernary's coverage of the
// analogous ConditionalOp path for the ConditionalStatement path.
TEST_P(SVUnitTests, MultiConditionIf)
{
  check_bmc("multi_cond_if.sv", 0, ProverResult::UNKNOWN);
}

// A `case` statement's `default:` arm only applies when no other item
// matched: process_statement's Case handler gives the default arm's
// condition as `condition AND NOT(any item matched)`, excluding every
// item's own match condition, rather than the bare outer `condition`.
TEST_P(SVUnitTests, CaseStatementDefaultOnlyWhenNoMatch)
{
  check_bmc("case_default.sv", 2);
}

// ---------------------------------------------------------------------------
// casex/casez wildcard matching: the `?` don't-care bits in `4'b1??1`
// make the case-item literal a 4-state value with unknown bits. The
// Case statement handler special-cases casex/casez: for a constant item
// pattern, it builds a (mask, value) pair from the pattern's own X
// (casex) or Z (both; `?` is an alias for `z`) bits and compares
// `(sel & mask) == value`, ignoring exactly the wildcard positions,
// instead of comparing the raw literal directly (which would otherwise
// reach expr_to_term()'s generic, wildcard-unaware IntegerLiteral case).
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, CasexCasezWildcard)
{
  check_bmc("casex_casez.sv", 4, ProverResult::UNKNOWN);
}

// `case (x) matches ... endcase` (StatementKind::PatternCase) is a
// distinct statement kind from plain `case`/`casex`/`casez`, matching
// each item's *pattern* rather than comparing values. A constant
// pattern is the degenerate case and pins the basic wiring.
TEST_P(SVUnitTests, PatternCase)
{
  check_prover<KInduction>("pattern_case.sv", 6, ProverResult::TRUE);
}

// A `.v` variable pattern binds what it matched for the arm to read,
// and matches anything -- so the items have to be first-match-wins
// (or the catch-all would overwrite the constant arm above it) and
// the case counts as exhaustive with no `default` present.
TEST_P(SVUnitTests, PatternCaseFirstMatchWins)
{
  check_prover<KInduction>(
      "pattern_case_first_match.sv", 6, ProverResult::TRUE);
}

// A structure pattern over a packed struct, mixing a constant field
// test with a binding one.
TEST_P(SVUnitTests, PatternCaseStructurePattern)
{
  check_prover<KInduction>("pattern_case_structure.sv", 6, ProverResult::TRUE);
}

// A `&&&` filter, which reads the name its own item's pattern bound
// and so needs that binding in scope while the guard is built.
TEST_P(SVUnitTests, PatternCaseFilter)
{
  check_prover<KInduction>("pattern_case_filter.sv", 6, ProverResult::TRUE);
}

// Not exhaustive and no `default`, so the target holds its value on
// the paths that assign nothing -- the latch synthesis infers.
TEST_P(SVUnitTests, PatternCaseIncompleteLatches)
{
  check_prover<KInduction>("pattern_case_incomplete.sv", 6, ProverResult::TRUE);
}

// A `tagged` pattern matches on a tagged union's discriminant, which
// a packed one keeps at a defined position (LRM 7.3.2): the tag in
// the top bits holding the member's declaration index, each member
// right-justified below.
TEST_P(SVUnitTests, PatternCaseTagged)
{
  check_prover<KInduction>("pattern_case_tagged.sv", 6, ProverResult::TRUE);
}

// The tag has to be what decides it. Reading the payload alone would
// make the last arm always win, which this refutes.
TEST_P(SVUnitTests, PatternCaseTaggedIgnoringTagFails)
{
  check_bmc("pattern_case_tagged_ignored.sv", 0, ProverResult::FALSE);
}

// ---------------------------------------------------------------------------
// `while`/`do-while`/`repeat`/`foreach` loop unrolling. A plain
// procedural scratch variable (`int i;` mutated by `i = i + 1;` inside
// the loop body) is neither a wire nor a state var, so writes to it go
// through slang's own constant evaluator and the loop_var_terms_ map
// instead of the normal wire/state-var write path. `while`/`do-while`
// conditions and `repeat` counts must be compile-time constants, same
// as `for` bounds (a runtime-dependent bound throws PonoException).
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, WhileLoop) { check_bmc("while_loop.sv", 2); }

TEST_P(SVUnitTests, DoWhileLoop) { check_bmc("do_while_loop.sv", 2); }

TEST_P(SVUnitTests, RepeatLoop) { check_bmc("repeat_loop.sv", 2); }

TEST_P(SVUnitTests, ForeachLoop) { check_bmc("foreach_loop.sv", 2); }

// ---------------------------------------------------------------------------
// `break`/`continue`/`disable`, scoped to compile-time-constant guard
// conditions: process_statement() throws a LoopControlSignal, caught
// by the nearest enclosing ForLoop (Break/Continue) or matching named
// Block (Disable). This relies on the `if` guarding each one being
// const-evaluable (each fixture below compares against an
// already-unrolled `for`-loop counter), so only the branch actually
// taken in C++ runs and the signal propagates correctly. See
// Gap_BreakRuntimeDependent below for the runtime-dependent case.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, BreakInForLoop)
{
  check_bmc("break_in_for.sv", 6, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, ContinueInForLoop)
{
  check_bmc("continue_in_for.sv", 6, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, DisableNamedBlock)
{
  check_bmc("disable_named_block.sv", 6, ProverResult::UNKNOWN);
}

// A `break` guarded by a runtime signal rather than a compile-time
// constant can't be modeled as C++-level control flow today (the
// general symbolic-guard path processes both `if` arms
// unconditionally) -- but a statically-bounded loop with a
// data-dependent early exit is ordinary synthesizable RTL (e.g. a
// priority-encoder-style search-with-break), not an inherent modeling
// impossibility like a truly unbounded `forever`. The fixture's own
// assertion demonstrates the correctly-encoded reachable behavior
// (`cond` staying low lets the loop run to completion and set
// `reg_any_set`); today this just throws instead.
TEST_P(SVUnitTests, Gap_BreakRuntimeDependent)
{
  check_bmc("break_runtime_dependent.sv", 2);
}

// Procedural immediate assertion (`assert (expr);`, distinct from
// `assert property (...)`): StatementKind::ImmediateAssertion builds a
// safety property (for `assert`) or a standing constraint (for
// `assume`/`restrict`), guarded by the accumulated path condition
// rather than treated as always-active.
//
// Unlike the register-based examples elsewhere in this suite, an
// immediate assertion inside always_ff checks *current-cycle* values with
// no register latency involved (rst and a are both plain inputs here, not
// registers) -- rst is already free (not forced) starting at cycle 1, so
// BMC can pick rst == 0 and a == 7 simultaneously at cycle 1 itself, not
// "one cycle after release" the way a registered value would need.
TEST_P(SVUnitTests, ImmediateAssert) { check_bmc("immediate_assert.sv", 1); }

// Immediate `assume` (Assume/Restrict share the same code path as assert).
TEST_P(SVUnitTests, ImmediateAssume)
{
  check_bmc("immediate_assume.sv", 4, ProverResult::UNKNOWN);
}

// A void task call used as a bare statement (`bump(a, b);`): the
// body is inlined and its output argument written back to the
// caller's variable.
TEST_P(SVUnitTests, VoidTaskCall)
{
  check_prover<KInduction>("void_task_call.sv", 8, ProverResult::TRUE);
}

// The argument directions the write-back has to distinguish, and a
// call on only one path, whose write-back must be guarded by that
// condition rather than applied unconditionally.
TEST_P(SVUnitTests, TaskCallArgumentDirections)
{
  check_prover<KInduction>("task_call_args.sv", 8, ProverResult::TRUE);
}

// What copy-in copy-out inlining cannot model.
TEST_P(SVUnitTests, TaskCallRecursiveRejected)
{
  expect_encode_throws("task_call_recursive.sv");
}

TEST_P(SVUnitTests, TaskCallRefArgRejected)
{
  expect_encode_throws("task_call_ref_arg.sv");
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVStatementsTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

// A variable whose only driver is an `initial` block: state rather
// than a free input, and given an update that holds, since a state
// var without one is treated as an input and would drift.

TEST_P(SVUnitTests, InitialOnlyDriver)
{
  check_prover<KInduction>("initial_only_driver.sv", 8, ProverResult::TRUE);
}

TEST_P(SVUnitTests, InitialOnlyDriverUnwrittenIsFree)
{
  check_bmc("initial_only_driver_fails.sv", 1);
}

// A runtime-indexed initial write needs no fixed slice: splicing
// onto the variable's own term pins the selected bits and leaves
// the rest of the initial value free.

TEST_P(SVUnitTests, InitialDynamicWrite)
{
  check_prover<KInduction>("initial_dynamic_write.sv", 8, ProverResult::TRUE);
}

TEST_P(SVUnitTests, InitialDynamicWriteRestIsFree)
{
  check_bmc("initial_dynamic_write_fails.sv", 0);
}

// A combinational block that assigns on some paths but not others
// infers a latch, as synthesis does -- per variable, so an
// A `case` covering every value of its selector with no `default`.
// The definite-assignment scan is syntactic and reads that as a path
// that writes nothing, so it marks the target a latch -- which would
// make a combinational signal a register and delay it a cycle.
// Whether the fallback to the old value is reachable is a question
// for the solver, and here it is not.
TEST_P(SVUnitTests, CaseFullCoverageIsCombinational)
{
  check_prover<KInduction>("case_full_coverage.sv", 6, ProverResult::TRUE);
}

// The same design, claiming what only the register reading would
// satisfy. Inferring a latch would *prove* this -- a false proof of
// behaviour the design does not have, which is the worse half of
// getting the previous test wrong.
TEST_P(SVUnitTests, CaseFullCoverageIsNotDelayed)
{
  check_bmc("case_full_coverage_delayed.sv", 1, ProverResult::FALSE);
}

// The shape it shows up as in practice: an enum state fully
// enumerated, no `default`.
TEST_P(SVUnitTests, CaseEnumFullCoverageIsCombinational)
{
  check_prover<KInduction>("case_enum_full_coverage.sv", 6, ProverResult::TRUE);
}

// unconditionally-assigned target in the same block stays a wire.

TEST_P(SVUnitTests, AlwaysCombLatch)
{
  check_prover<KInduction>("always_comb_latch.sv", 12, ProverResult::TRUE);
}

TEST_P(SVUnitTests, AlwaysCombLatchIsNotAWire)
{
  check_bmc("always_comb_latch_fails.sv", 1);
}

TEST_P(SVUnitTests, AlwaysCombFullyAssigned)
{
  check_prover<KInduction>("always_comb_full.sv", 12, ProverResult::TRUE);
}

}  // namespace pono_tests
