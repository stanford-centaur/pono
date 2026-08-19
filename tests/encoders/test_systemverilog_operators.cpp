#ifdef WITH_SLANG

#include "sv_test_fixture.h"

using namespace pono;
using namespace smt;

namespace pono_tests {

// ---------------------------------------------------------------------------
// Arithmetic + comparison, bitwise + logical, shift + unary + reduction,
// ternary + select + concatenation: each pair below replaces what used to
// be several near-duplicate bare single-operator tests with one algebraic
// invariant checked against free inputs across many cycles (a real proof,
// not a single hand-picked trace) plus a companion that deliberately
// breaks the invariant.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, ArithCompareHolds)
{
  check_bmc("arith_compare.sv", 8, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, ArithCompareFails)
{
  check_bmc("arith_compare_fails.sv", 1);
}

TEST_P(SVUnitTests, BitwiseLogicalHolds)
{
  check_bmc("bitwise_logical.sv", 6, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, BitwiseLogicalFails)
{
  check_bmc("bitwise_logical_fails.sv", 1);
}

TEST_P(SVUnitTests, ShiftUnaryReductionHolds)
{
  check_bmc("shift_unary_reduction.sv", 6, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, ShiftUnaryReductionFails)
{
  check_bmc("shift_unary_reduction_fails.sv", 1);
}

TEST_P(SVUnitTests, TernarySelectConcatHolds)
{
  check_bmc("ternary_select_concat.sv", 6, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, TernarySelectConcatFails)
{
  check_bmc("ternary_select_concat_fails.sv", 1);
}

// ---------------------------------------------------------------------------
// Replication + every sized-literal spelling; the three equivalent
// combinational-block styles (assign / always_comb / legacy always @*).
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, ReplicationSizedLiteralsHold)
{
  check_bmc("replication_sized_literals.sv", 4, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, CombPathsEquivalenceHolds)
{
  check_bmc("comb_paths.sv", 4, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, CombPathsEquivalenceFails)
{
  check_bmc("comb_paths_fails.sv", 2);
}

// ---------------------------------------------------------------------------
// Operators absent from the encoder's BinaryOperator/UnaryOperator
// switches -- confirmed by inspection to hit the `throw PonoException`
// default case in expr_to_term(), so each of these either passes as a
// correctness check (if support exists) or fails via an uncaught
// exception (if not) -- verified per-case, not assumed.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, Gap_PowerOperator) { check_bmc("power_op.sv", 2); }

TEST_P(SVUnitTests, Gap_EqualityVariants)
{
  check_bmc("equality_variants.sv", 4, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, Gap_ReductionNandNorXnor)
{
  check_bmc("reduction_nand_nor_xnor.sv", 4, ProverResult::UNKNOWN);
}

// FIXED: process_statement()'s ExpressionStatement case only ever
// recognized ExpressionKind::Assignment, so a standalone `i++;`/`--i;`
// statement (distinct from the same operators used as a `for`-loop
// step, which slang's own constant evaluator already handled via a
// completely separate path) fell through and did nothing. Added a
// UnaryOp branch that reuses the same lvalue-resolution/commit
// machinery as plain assignment (refactored into begin_write()/
// commit_write() local lambdas shared by both), reading the current
// value directly instead of evaluating an RHS expression -- so this
// works for any lvalue shape resolve_lvalue() supports, not just a
// plain scalar (see ElementIncrement/StructFieldIncrement below).
TEST_P(SVUnitTests, IncrementDecrement) { check_bmc("inc_dec.sv", 4); }

// `++` on a packed-array element and a packed-struct field,
// confirming the fix above is genuinely general rather than a
// scalar-only special case.
TEST_P(SVUnitTests, ElementIncrement) { check_bmc("element_increment.sv", 4); }

TEST_P(SVUnitTests, StructFieldIncrement)
{
  check_bmc("struct_field_increment.sv", 4);
}

TEST_P(SVUnitTests, Gap_StreamingOperator) { check_bmc("streaming_op.sv", 2); }

TEST_P(SVUnitTests, Gap_UnaryPlusIdentity) { check_bmc("unary_plus.sv", 2); }

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVOperatorsTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests

#endif
