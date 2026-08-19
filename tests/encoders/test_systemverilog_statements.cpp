#ifdef WITH_SLANG

#include "sv_test_fixture.h"

using namespace pono;
using namespace smt;

namespace pono_tests {

// ---------------------------------------------------------------------------
// Migrated from the original test_systemverilog.cpp -- already compositional
// (procedural for-loop + compound assignment / element-select LHS composed
// with always_ff), kept as-is per the triage in the plan.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, CompoundAssignOrReduce)
{
  check_bmc("compound_assign.sv", 2);
}

TEST_P(SVUnitTests, ElementSelectLhs) { check_bmc("element_select_lhs.sv", 2); }

// Basic block-level smoke checks -- module ports + always_ff (counter.sv)
// and `initial` blocks pinning initial state (initial_block.sv) -- kept
// from the original suite verbatim rather than folded away: the pattern
// they exercise (a bare always_ff counter; a design with no `rst` port at
// all, relying solely on `initial` to pin state) recurs as a building
// block throughout this whole suite, so it's worth keeping one minimal,
// direct test of each rather than only ever seeing it composed with
// something else.
TEST_P(SVUnitTests, EncodeCounter) { check_bmc("counter.sv", 5); }

TEST_P(SVUnitTests, InitialBlockSetsState) { check_bmc("initial_block.sv", 0); }

// ---------------------------------------------------------------------------
// `priority if` / `unique case` as semantic modifiers, replacing the old
// bare if_else.sv / case_stmt.sv with one denser, paired test.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, PriorityIfUniqueCaseFails)
{
  check_bmc("priority_if_unique_case.sv", 2);
}

TEST_P(SVUnitTests, PriorityIfUniqueCaseHolds)
{
  check_bmc("priority_if_unique_case_holds.sv", 6, ProverResult::UNKNOWN);
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

// ---------------------------------------------------------------------------
// `while`/`do-while`/`repeat`/`foreach` loop unrolling. A plain
// procedural scratch variable (`int i;` mutated by `i = i + 1;` inside
// the loop body) needs its own fast path in ExpressionStatement's
// Assignment/++/-- handling: a write to any symbol currently bound as a
// compile-time-unrolled local is evaluated by slang's own constant
// evaluator and mirrored back into the same loop_var_terms_ map
// `for`-loop counters use, rather than the wire/state-var write path
// (which only applies to real design signals). `while`/`do-while`
// conditions and `repeat` counts follow the same compile-time-constant-
// only contract as `for` bounds and break/continue/disable (clear
// PonoException, not a wrong verdict, for a runtime-dependent bound);
// `foreach` is scoped to a single statically-sized dimension with a
// loop variable, which is what `foreach (arr[i])` over a fixed-size
// packed vector produces.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, WhileLoop) { check_bmc("while_loop.sv", 2); }

TEST_P(SVUnitTests, DoWhileLoop) { check_bmc("do_while_loop.sv", 2); }

TEST_P(SVUnitTests, RepeatLoop) { check_bmc("repeat_loop.sv", 2); }

TEST_P(SVUnitTests, ForeachLoop) { check_bmc("foreach_loop.sv", 2); }

// ---------------------------------------------------------------------------
// `break`/`continue`/`disable`, scoped to compile-time-constant
// conditions: process_statement() throws a LoopControlSignal from
// Break/Continue/Disable, caught by the nearest enclosing ForLoop
// (Break/Continue) or matching named Block (Disable). This relies on
// the Conditional case const-evaluating an `if`'s condition first
// (falling back to the existing symbolic-guard path when that fails)
// -- each fixture below guards its break/continue/disable with a
// comparison against an already-unrolled `for`-loop counter, so the
// const-eval succeeds and the signal can propagate out of the one
// branch actually taken in C++, the same way real control flow would.
// See Unsupported_BreakRuntimeDependent below for the (deliberately
// unsupported, clean-error) runtime-dependent case.
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
// constant can't be modeled as C++-level control flow at all (the
// general symbolic-guard path processes both `if` arms
// unconditionally) -- a genuine architectural boundary of the
// compile-time-unrolling model, not a "not implemented yet" gap, so
// it's named like the other deliberate-non-goal Unsupported_ cases.
// Must throw a clear PonoException rather than silently doing nothing
// or applying the wrong verdict either way.
TEST_P(SVUnitTests, Unsupported_BreakRuntimeDependent)
{
  expect_encode_throws("break_runtime_dependent.sv");
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

// Immediate `assume` (Assume/Restrict share the same fixed code path).
TEST_P(SVUnitTests, ImmediateAssume)
{
  check_bmc("immediate_assume.sv", 4, ProverResult::UNKNOWN);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVStatementsTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests

#endif
