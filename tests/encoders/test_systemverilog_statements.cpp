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

// FIXED: a `case` statement's `default:` arm used to apply
// unconditionally alongside whichever other item already matched
// (process_statement's Case handler processed the default arm with
// the bare outer `condition` rather than `condition AND NOT(any item
// matched)`), so any case statement with a `default:` clause always
// ended up applying the default's assignment regardless of which item
// actually matched. Now the default arm's condition explicitly
// excludes every item's match condition.
TEST_P(SVUnitTests, CaseStatementDefaultOnlyWhenNoMatch)
{
  check_bmc("case_default.sv", 2);
}

// ---------------------------------------------------------------------------
// FIXED: the `?` don't-care bits in `4'b1??1` make the case-item literal a
// 4-state value with unknown bits, which used to get handed straight to
// expr_to_term()'s generic (wildcard-unaware) IntegerLiteral case, which in
// turn handed an X-containing decimal string straight to the solver and
// crashed with a raw solver-library exception ("invalid decimal string" /
// "mpz_set_str") rather than implementing (or even cleanly rejecting)
// wildcard matching. The Case statement handler now special-cases
// casex/casez: for a constant item pattern, it builds a (mask, value) pair
// from the pattern's own X (casex) or Z (both; `?` is an alias for `z`)
// bits and compares `(sel & mask) == value`, ignoring exactly the wildcard
// positions, instead of comparing the raw (invalid) literal directly.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, CasexCasezWildcard)
{
  check_bmc("casex_casez.sv", 4, ProverResult::UNKNOWN);
}

// ---------------------------------------------------------------------------
// Loop kinds beyond the already-supported `for`.  process_statement()'s
// default case (confirmed by inspection) silently skips any statement kind
// it doesn't recognize, logging only a level-1 warning -- so an
// unimplemented loop kind doesn't throw, it just makes the loop body a
// no-op.  Each of these mirrors the OR-reduction compound_assign.sv/
// while_loop.sv already prove works via `for`, so a Gap_ result here
// isolates the loop *keyword* as the point of failure.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, Gap_WhileLoop) { check_bmc("while_loop.sv", 2); }

TEST_P(SVUnitTests, Gap_DoWhileLoop) { check_bmc("do_while_loop.sv", 2); }

TEST_P(SVUnitTests, Gap_RepeatLoop) { check_bmc("repeat_loop.sv", 2); }

TEST_P(SVUnitTests, Gap_ForeachLoop) { check_bmc("foreach_loop.sv", 2); }

// ---------------------------------------------------------------------------
// FIXED, scoped to compile-time-constant conditions: process_statement()
// gained a LoopControlSignal thrown by Break/Continue/Disable and caught by
// the nearest enclosing ForLoop (Break/Continue) or matching named Block
// (Disable). This only works because the Conditional case now also tries
// to const-evaluate an `if`'s condition first (falling back to the
// existing symbolic-guard path when that fails) -- each fixture below
// guards its break/continue/disable with a comparison against an
// already-unrolled `for`-loop counter, so the const-eval succeeds and the
// signal can propagate out of the one branch actually taken in C++, the
// same way real control flow would.  See Gap_BreakRuntimeDependent below
// for the (deliberately unsupported, clean-error) runtime-dependent case.
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

// FIXED: procedural immediate assertion (`assert (expr);`, distinct from
// `assert property (...)`) had no StatementKind::ImmediateAssertion case at
// all and so produced zero properties, contradicting the encoder header's
// own doc-comment claim that "Basic SVA immediate assertions" are
// supported. Added a case that builds a safety property (for `assert`) or
// a standing constraint (for `assume`/`restrict`), guarded by the
// accumulated path condition rather than treated as always-active.
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
