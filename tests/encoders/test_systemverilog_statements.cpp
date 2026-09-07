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

// A constant element-select lvalue whose index is out of range for its
// base (`flag[10]` into a 4-bit `flag`) -- the LRM permits this
// (writes are a no-op, reads return 'x), but this encoder has no such
// semantics.
TEST_P(SVUnitTests, Gap_ElementSelectOutOfBoundsLhs)
{
  check_bmc("element_select_out_of_bounds_lhs.sv", 2, ProverResult::UNKNOWN);
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
  SystemVerilogEncoder::encode(sv_path("concat_lhs.sv"), fts);
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
// (`{>>{hi, lo}} <= a;`) is ExpressionKind::Streaming, distinct from a
// plain concatenation-target LHS (ExpressionKind::Concatenation,
// already supported above). resolve_lvalue() has no case for it at
// all.
TEST_P(SVUnitTests, Gap_StreamingConcatLhs)
{
  check_bmc("streaming_concat_lhs.sv", 2, ProverResult::UNKNOWN);
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
// distinct statement kind from plain case/casex/casez
// (StatementKind::Case) that pre_scan_state_vars()'s
// collect_blocking_targets()/collect_nonblocking_targets() don't
// recognize either -- but since process_statement() itself also
// doesn't process it (falling to the generic unhandled-statement-kind
// default), the two omissions are consistent: no write inside it is
// ever pre-scanned *or* applied, and the skip is logged. A real
// mainstream-RTL gap, not a deliberate non-goal.
TEST_P(SVUnitTests, Gap_PatternCase)
{
  check_bmc("pattern_case.sv", 2, ProverResult::UNKNOWN);
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

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVStatementsTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
