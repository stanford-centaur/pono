#include "engines/kinduction.h"
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

// A comparison is encoded in whichever sort its context wants: a
// native Bool where only its truth value is needed, a 1-bit bit-vector
// where an actual value is (assigned into a net, fed to a concat).
// These two check that the two forms never disagree -- the holds
// variant uses each comparison in both roles at once.
TEST_P(SVUnitTests, PredicateValueAndConditionHolds)
{
  check_bmc("predicate_value_and_condition.sv", 6, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, PredicateValueAndConditionFails)
{
  check_bmc("predicate_value_and_condition_fails.sv", 1);
}

// A $past call inside a `&&` operand must be converted exactly once:
// its history chain is not memoized, so a second conversion would add
// a second chain of latches tracking the same value.  That is
// invisible to any property verdict, so check the latch count.
TEST_P(SVUnitTests, PastInLogicalAndBuildsOneChain)
{
  SmtSolver s = create_solver(GetParam());
  FunctionalTransitionSystem fts(s);
  SystemVerilogEncoder::encode(fts, sv_path("past_in_logical_and.sv"));

  size_t chain_latches = 0;
  for (const auto & sv : fts.statevars()) {
    if (sv->to_string().find("__sva_past_") != std::string::npos) {
      ++chain_latches;
    }
  }
  // One `$past(a)`, one cycle of delay, so exactly one latch.
  EXPECT_EQ(chain_latches, 1u);
}

TEST_P(SVUnitTests, PastInLogicalAndHolds)
{
  check_bmc("past_in_logical_and.sv", 6, ProverResult::UNKNOWN);
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
// Power ('**'), case equality ('===='/'!=='), wildcard equality ('==?'/
// '!=?'), and reduction NAND/NOR/XNOR ('~&'/'~|'/'~^').
//   - Power: scoped to a compile-time-constant exponent, unrolled into
//     repeated multiplication -- a non-constant exponent throws (real BV
//     exponentiation isn't part of the SMT BV theory).
//   - Case equality: identical to logical equality/inequality, since this
//     encoder's pure-BV model has no X/Z to make them actually differ.
//   - Wildcard equality: the same (mask, value) technique as casex/casez,
//     applied to the right operand's X/Z bits per the LRM; falls back to
//     plain equality for a non-constant right operand.
//   - Reduction NAND/NOR/XNOR: the existing AND/OR/XOR reduction logic,
//     negated.
// equality_variants.sv's wildcard-equality case has an X-containing literal
// on the right operand, which is special-cased before the encoder's eager
// BinaryOp operand conversion so it never reaches the generic (wildcard-
// unaware) integer-literal path.
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

// `i++;`/`--i;` as a standalone statement (distinct from the same
// operator used as a `for`-loop step, which slang's constant evaluator
// handles separately). Reuses the assignment lvalue-resolution/commit
// machinery, so this works for any lvalue shape resolve_lvalue()
// supports, not just a plain scalar (see ElementIncrement/
// StructFieldIncrement below).
TEST_P(SVUnitTests, IncrementDecrement) { check_bmc("inc_dec.sv", 4); }

// `++` on a packed-array element and a packed-struct field, confirming
// the general lvalue support above rather than a scalar-only special
// case.
TEST_P(SVUnitTests, ElementIncrement) { check_bmc("element_increment.sv", 4); }

TEST_P(SVUnitTests, StructFieldIncrement)
{
  check_bmc("struct_field_increment.sv", 4);
}

// Streaming concatenation ('{<<{...}}'/'{>>{...}}'), scoped to a single
// stream with no `with` sub-range (real usage -- reversing/regrouping
// one packed value's bits or byte-lanes), reassembling slice-sized
// chunks in reverse order via Extract+Concat; the LRM's full generality
// (multiple streams, `with` ranges, dynamically-sized queues) throws a
// clear error.
TEST_P(SVUnitTests, StreamingOperator) { check_bmc("streaming_op.sv", 2); }

// Unary '+' (a no-op per the LRM).
TEST_P(SVUnitTests, UnaryPlusIdentity) { check_bmc("unary_plus.sv", 2); }

// A `&&&`-joined multi-condition ternary (`a &&& b ? x : y`, LRM
// 11.4.11) is legal outside case/if context too; expr_to_term() must
// AND together every joined condition, not just the first.
TEST_P(SVUnitTests, MultiConditionTernary)
{
  check_bmc("multi_cond_ternary.sv", 0, ProverResult::UNKNOWN);
}

// A user-defined function called with a symbolic argument: the body
// is inlined where the call appears. Proved rather than unrefuted,
// so the value it produces is pinned, not just its encodability.
TEST_P(SVUnitTests, UserFunctionCall)
{
  check_prover<KInduction>("user_function_call.sv", 4, ProverResult::TRUE);
}

// A body that is more than one return: locals, assignment to the
// function's own name, a conditional overwrite, an unrolled loop.
TEST_P(SVUnitTests, UserFunctionBody)
{
  check_prover<KInduction>("user_function_body.sv", 4, ProverResult::TRUE);
}

// Repeated, nested, and function-to-function calls, each of which
// rebinds the same formals.
TEST_P(SVUnitTests, UserFunctionNested)
{
  check_prover<KInduction>("user_function_nested.sv", 4, ProverResult::TRUE);
}

// What inlining cannot reach: a call that never stops expanding, a
// return whose value depends on which path ran, and a call that
// writes back to its caller.
TEST_P(SVUnitTests, UserFunctionRecursiveRejected)
{
  expect_encode_throws("user_function_recursive.sv");
}

TEST_P(SVUnitTests, UserFunctionCondReturnRejected)
{
  expect_encode_throws("user_function_cond_return.sv");
}

TEST_P(SVUnitTests, UserFunctionOutputArgRejected)
{
  expect_encode_throws("user_function_output_arg.sv");
}

// A wildcard pattern wider than 64 bits, via `==?`, `!=?` and a
// casex arm. Each is asserted equivalent to comparing exactly the
// bits the pattern pins; the low- and high-byte variants together
// pin the bit order.
TEST_P(SVUnitTests, WildcardWidePattern)
{
  check_prover<KInduction>("wildcard_wide_pattern.sv", 4, ProverResult::TRUE);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVOperatorsTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

// ---------------------------------------------------------------------------
// `x` and `z` bits in a literal. Modelled as unconstrained bits, in
// keeping with how this encoder treats X everywhere else, rather than
// reaching the solver as a decimal string containing an `x`.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, XzLiteralKnownBits)
{
  check_prover<KInduction>("xz_literal.sv", 12, ProverResult::TRUE);
}

TEST_P(SVUnitTests, XzLiteralUnknownBitsAreFree)
{
  check_bmc("xz_literal_fails.sv", 1);
}

TEST_P(SVUnitTests, XzLiteralsAreIndependent)
{
  check_bmc("xz_literal_independent.sv", 1);
}

}  // namespace pono_tests
