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

// ---------------------------------------------------------------------------
// GAP (confirmed empirically): the `?` don't-care bits in `4'b1??1` make
// the case-item literal a 4-state value with unknown bits.  The encoder
// doesn't special-case that when converting the literal to an SMT term, so
// it hands an X-containing decimal string straight to the solver and
// crashes with a raw solver-library exception ("invalid decimal string" /
// "mpz_set_str") rather than a clean PonoException -- itself worth noting
// alongside the missing wildcard-matching semantics.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, Gap_CasexCasezWildcard)
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
// GAP (confirmed empirically): break / continue / disable inside the
// already-supported `for` loop are all silently skipped rather than
// honored -- each fixture's distinguishing construction (see the .sv
// comments) observed the "control-flow statement ignored" outcome
// (property falsified) rather than the "honored" outcome (property holds)
// asserted here per correct LRM semantics.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, Gap_BreakInForLoop)
{
  check_bmc("break_in_for.sv", 6, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, Gap_ContinueInForLoop)
{
  check_bmc("continue_in_for.sv", 6, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, Gap_DisableNamedBlock)
{
  check_bmc("disable_named_block.sv", 6, ProverResult::UNKNOWN);
}

// ---------------------------------------------------------------------------
// GAP (confirmed empirically, and contradicting the encoder header's own
// doc-comment claim that "Basic SVA immediate assertions" are supported):
// procedural immediate assertion (`assert (expr);`, distinct from
// `assert property (...)`) produces zero properties -- enc.propvec() is
// empty, so the ASSERT_EQ(propvec().size(), 1u) inside check_bmc() fails
// before BMC even runs.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, Gap_ImmediateAssert)
{
  check_bmc("immediate_assert.sv", 2);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVStatementsTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests

#endif
