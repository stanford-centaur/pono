#include "engines/kinduction.h"
#include "sv_test_fixture.h"

using namespace pono;
using namespace smt;

namespace pono_tests {

// ---------------------------------------------------------------------------
// Expressions and operators.
//
// Most operator families are covered by an algebraic invariant checked
// against free inputs over many cycles, rather than by a hand-picked
// trace, and each such invariant is paired with a design that breaks it.
// The pair is the point: the holding half alone cannot tell a correct
// encoding from one that makes the property unreachable.
//
// Grouped by the part of expression handling under test. A test for a
// new operator belongs with the family it joins; one for how an
// expression is *placed* (its sort, its lvalue, its literal form)
// belongs in the matching section further down.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Operator families, one invariant plus one broken twin each.
// ---------------------------------------------------------------------------

// Arithmetic and comparison.
TEST_P(SVUnitTests, ArithCompareHolds)
{
  check_bmc("arith_compare.sv", 8, ProverResult::UNKNOWN);
}

// The same invariant perturbed, so it must be refuted.
TEST_P(SVUnitTests, ArithCompareFails)
{
  check_bmc("arith_compare_fails.sv", 1);
}

// Bitwise and logical.
TEST_P(SVUnitTests, BitwiseLogicalHolds)
{
  check_bmc("bitwise_logical.sv", 6, ProverResult::UNKNOWN);
}

// The same invariant perturbed.
TEST_P(SVUnitTests, BitwiseLogicalFails)
{
  check_bmc("bitwise_logical_fails.sv", 1);
}

// Shift, unary and reduction.
TEST_P(SVUnitTests, ShiftUnaryReductionHolds)
{
  check_bmc("shift_unary_reduction.sv", 6, ProverResult::UNKNOWN);
}

// The same invariant perturbed.
TEST_P(SVUnitTests, ShiftUnaryReductionFails)
{
  check_bmc("shift_unary_reduction_fails.sv", 1);
}

// Ternary, select and concatenation.
TEST_P(SVUnitTests, TernarySelectConcatHolds)
{
  check_bmc("ternary_select_concat.sv", 6, ProverResult::UNKNOWN);
}

// The same invariant perturbed.
TEST_P(SVUnitTests, TernarySelectConcatFails)
{
  check_bmc("ternary_select_concat_fails.sv", 1);
}

// Replication together with every sized-literal spelling.
TEST_P(SVUnitTests, ReplicationSizedLiteralsHold)
{
  check_bmc("replication_sized_literals.sv", 4, ProverResult::UNKNOWN);
}

// The three equivalent combinational-block styles -- `assign`,
// `always_comb`, and legacy `always @*` -- computing the same function.
TEST_P(SVUnitTests, CombPathsEquivalenceHolds)
{
  check_bmc("comb_paths.sv", 4, ProverResult::UNKNOWN);
}

// One style computing something else, so the three disagree.
TEST_P(SVUnitTests, CombPathsEquivalenceFails)
{
  check_bmc("comb_paths_fails.sv", 2);
}

// Reduction NAND, NOR and XNOR (`~&`, `~|`, `~^`).
TEST_P(SVUnitTests, ReductionNandNorXnor)
{
  check_bmc("reduction_nand_nor_xnor.sv", 4, ProverResult::UNKNOWN);
}

// Power (`**`) with a compile-time-constant exponent, unrolled into
// repeated multiplication. A non-constant exponent is refused, since
// the SMT bit-vector theory has no exponentiation.
TEST_P(SVUnitTests, PowerOperator) { check_bmc("power_op.sv", 2); }

// Unary `+`, a no-op per the LRM.
TEST_P(SVUnitTests, UnaryPlusIdentity) { check_bmc("unary_plus.sv", 2); }

// A streaming concatenation (`{<<{...}}` / `{>>{...}}`) as an
// expression, reassembling slice-sized chunks in reverse order.
TEST_P(SVUnitTests, StreamingOperator) { check_bmc("streaming_op.sv", 2); }

// ---------------------------------------------------------------------------
// Equality forms. Case equality (`===`/`!==`) coincides with logical
// equality here, since a pure bit-vector model has no X or Z for them
// to differ over; wildcard equality (`==?`/`!=?`) masks the right
// operand's X and Z bits, the same rule casex and casez use.
// ---------------------------------------------------------------------------

// All four spellings together, with an X-containing literal on the
// right of the wildcard comparison.
TEST_P(SVUnitTests, EqualityVariants)
{
  check_bmc("equality_variants.sv", 4, ProverResult::UNKNOWN);
}

// A wildcard pattern wider than 64 bits, through `==?`, `!=?` and a
// casex arm. Each is checked against a comparison of exactly the bits
// the pattern pins, and the low- and high-byte variants together pin
// the bit order.
TEST_P(SVUnitTests, WildcardWidePattern)
{
  check_prover<KInduction>("wildcard_wide_pattern.sv", 4, ProverResult::TRUE);
}

// ---------------------------------------------------------------------------
// Set membership, `inside` (LRM 11.4.13): an or-reduction over the
// members, where an integral member is matched by wildcard equality and
// a `[lo:hi]` member by an inclusive range. The expression form and
// `case ... inside` share the rule, so both are checked against it.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, InsideExpression)
{
  check_prover<KInduction>("inside_expression.sv", 8, ProverResult::TRUE);
}

// A value outside the set.
TEST_P(SVUnitTests, InsideExpressionFails)
{
  check_bmc("inside_expression_fails.sv", 0);
}

// The `case ... inside` form, including a `[lo:hi]` member.
TEST_P(SVUnitTests, CaseInsideRange)
{
  check_prover<KInduction>("case_inside_range.sv", 12, ProverResult::TRUE);
}

// A selector outside every arm's set.
TEST_P(SVUnitTests, CaseInsideRangeFails)
{
  check_bmc("case_inside_range_fails.sv", 2);
}

// ---------------------------------------------------------------------------
// How an expression is placed: the sort its context asks for, and the
// conditions a ternary can carry.
// ---------------------------------------------------------------------------

// A comparison is encoded in whichever sort its context wants -- a
// native Bool where only its truth value is needed, a 1-bit bit-vector
// where a value is (assigned to a net, fed to a concatenation). These
// check the two forms agree, with the holding half using each
// comparison in both roles at once.
TEST_P(SVUnitTests, PredicateValueAndConditionHolds)
{
  check_bmc("predicate_value_and_condition.sv", 6, ProverResult::UNKNOWN);
}

// The two roles made to disagree.
TEST_P(SVUnitTests, PredicateValueAndConditionFails)
{
  check_bmc("predicate_value_and_condition_fails.sv", 1);
}

// A `&&&`-joined multi-condition ternary (`a &&& b ? x : y`, LRM
// 11.4.11), which is legal outside case and if context: every joined
// condition is ANDed together, not just the first.
TEST_P(SVUnitTests, MultiConditionTernary)
{
  check_bmc("multi_cond_ternary.sv", 0, ProverResult::UNKNOWN);
}

// ---------------------------------------------------------------------------
// `$past` inside an expression. Its history chain is not memoized, so
// how many times the call is converted is observable in the state
// vars even where no property verdict changes -- which is what the
// first test here looks at directly.
// ---------------------------------------------------------------------------

// One `$past(a)` with one cycle of delay builds exactly one latch,
// however many operands of the enclosing `&&` reach it.
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
  EXPECT_EQ(chain_latches, 1u);
}

// The same design's own invariant.
TEST_P(SVUnitTests, PastInLogicalAndHolds)
{
  check_bmc("past_in_logical_and.sv", 6, ProverResult::UNKNOWN);
}

// ---------------------------------------------------------------------------
// Increment and decrement, which reach the assignment lvalue machinery
// rather than being a scalar-only special case -- so the shapes below
// widen from a plain variable outwards.
// ---------------------------------------------------------------------------

// `i++` / `--i` as a standalone statement, distinct from the same
// operator as a `for`-loop step, which slang's constant evaluator
// handles separately.
TEST_P(SVUnitTests, IncrementDecrement) { check_bmc("inc_dec.sv", 4); }

// `++` on a packed-array element.
TEST_P(SVUnitTests, ElementIncrement) { check_bmc("element_increment.sv", 4); }

// `++` on a packed-struct field.
TEST_P(SVUnitTests, StructFieldIncrement)
{
  check_bmc("struct_field_increment.sv", 4);
}

// ---------------------------------------------------------------------------
// `x` and `z` bits in a literal, modelled as unconstrained bits in
// keeping with how this encoder treats X everywhere else.
// ---------------------------------------------------------------------------

// The known bits of such a literal keep their values.
TEST_P(SVUnitTests, XzLiteralKnownBits)
{
  check_prover<KInduction>("xz_literal.sv", 12, ProverResult::TRUE);
}

// The unknown bits take either value.
TEST_P(SVUnitTests, XzLiteralUnknownBitsAreFree)
{
  check_bmc("xz_literal_fails.sv", 1);
}

// Two X bits are independent of each other rather than one shared
// value.
TEST_P(SVUnitTests, XzLiteralsAreIndependent)
{
  check_bmc("xz_literal_independent.sv", 1);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVOperatorsTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
