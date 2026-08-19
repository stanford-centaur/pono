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
// FIXED: these operators were absent from the encoder's BinaryOperator/
// UnaryOperator switches, hitting the generic `throw PonoException` default
// case in expr_to_term().
//   - Power ('**'): scoped to a compile-time-constant exponent, unrolled
//     into repeated multiplication -- a non-constant exponent still throws
//     (real BV exponentiation isn't part of the SMT BV theory).
//   - Case equality ('===='/'!=='): identical to logical equality/
//     inequality, since this encoder's pure-BV model has no X/Z to make
//     them actually differ.
//   - Wildcard equality ('==?'/'!=?'): the same (mask, value) technique as
//     casex/casez, applied to the right operand's X/Z bits per the LRM;
//     falls back to plain equality for a non-constant right operand.
//   - Reduction NAND/NOR/XNOR ('~&'/'~|'/'~^'): the existing AND/OR/XOR
//     reduction logic, negated.
// equality_variants.sv's own draft had briefly hit the same "X-containing
// literal reaches the generic IntegerLiteral case" crash casex/casez had,
// since the encoder eagerly converts *both* BinaryOp operands before
// dispatching on the operator -- WildcardEquality/Inequality are now
// special-cased before that eager conversion.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, PowerOperator) { check_bmc("power_op.sv", 2); }

TEST_P(SVUnitTests, EqualityVariants)
{
  check_bmc("equality_variants.sv", 4, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, ReductionNandNorXnor)
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

// FIXED: streaming concatenation ('{<<{...}}'/'{>>{...}}') was a
// completely unhandled ExpressionKind. Scoped to a single stream with no
// `with` sub-range (real usage -- reversing/regrouping one packed value's
// bits or byte-lanes), reassembling slice-sized chunks in reverse order via
// Extract+Concat; the LRM's full generality (multiple streams, `with`
// ranges, dynamically-sized queues) still throws a clear error.
TEST_P(SVUnitTests, StreamingOperator) { check_bmc("streaming_op.sv", 2); }

// FIXED: unary '+' (a no-op per the LRM) was absent from the
// UnaryOperator switch.
TEST_P(SVUnitTests, UnaryPlusIdentity) { check_bmc("unary_plus.sv", 2); }

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVOperatorsTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests

#endif
