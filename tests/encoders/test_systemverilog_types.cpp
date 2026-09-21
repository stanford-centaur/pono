#include "engines/kinduction.h"
#include "sv_test_fixture.h"

using namespace pono;
using namespace smt;

namespace pono_tests {

// ---------------------------------------------------------------------------
// Packed multi-dimensional arrays (array-of-vectors): constant- and
// variable-index reads/writes on a packed array.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, PackedArrayConstIndexChain)
{
  check_bmc("array_const_index.sv", 7);
}

TEST_P(SVUnitTests, PackedArrayDynIndexReadWrite)
{
  check_bmc("array_dyn_index.sv", 16);
}

// Declared ranges other than `[n:0]`: ascending, not reaching zero,
// and crossing zero. The index is not the bit offset in any of them.

TEST_P(SVUnitTests, PackedArrayRanges)
{
  check_prover<KInduction>("packed_array_ranges.sv", 12, ProverResult::TRUE);
}

TEST_P(SVUnitTests, PackedArrayRangesFails)
{
  check_bmc("packed_array_ranges_fails.sv", 2);
}

// ---------------------------------------------------------------------------
// Unpacked arrays (register files / small memories), distinct from the
// packed-array tests above: these become a genuine SMT array sort, so
// element access is Select/Store rather than bit arithmetic.
//
// Supported as registers internal to one module; the exclusions below
// are asserted rather than assumed.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, UnpackedRegfileMemory)
{
  check_bmc("unpacked_regfile.sv", 3, ProverResult::UNKNOWN);
}

// The read-after-write invariant is 1-inductive, so BMC alone can only
// fail to refute it. Proving it is what distinguishes a correct
// encoding from an over-constrained one that happens to look quiet.
TEST_P(SVUnitTests, UnpackedRegfileMemoryProvable)
{
  check_prover<KInduction>("unpacked_regfile.sv", 12, ProverResult::TRUE);
}

// The other half of that guard: a memory does not preserve the address
// *next* to the one written, so this must be refuted. A vacuous or
// over-constrained array encoding would report no violation.
TEST_P(SVUnitTests, UnpackedRegfileMemoryFails)
{
  check_bmc("unpacked_regfile_fails.sv", 2);
}

// `foreach` reset, whole-array constant assignment, and a non-zero-based
// declared range (`mem[3:18]`, where element 0 of the SMT array is
// mem[3]) in one design.
TEST_P(SVUnitTests, UnpackedArrayReset)
{
  check_prover<KInduction>("unpacked_array_reset.sv", 12, ProverResult::TRUE);
}

// A bit select, a part select and a packed-struct field write, each
// inside one array element. Proving the read-back is what rules out
// the write having been dropped, which for an array would leave it
// free rather than merely stale.
TEST_P(SVUnitTests, UnpackedArraySubelementWrite)
{
  check_prover<KInduction>(
      "unpacked_array_subelement.sv", 12, ProverResult::TRUE);
}

// Ite over two array terms, which has no bit width to match.
TEST_P(SVUnitTests, UnpackedArraySelect)
{
  check_prover<KInduction>("unpacked_array_select.sv", 12, ProverResult::TRUE);
}

// An index outside the declared range: the write is ignored and the
// read is X, and neither may touch a real cell. The paired
// refutation below is what rules out the phantom cell the truncated
// index used to create.
TEST_P(SVUnitTests, UnpackedArrayOutOfRange)
{
  check_prover<KInduction>(
      "unpacked_array_out_of_range.sv", 12, ProverResult::TRUE);
}

TEST_P(SVUnitTests, UnpackedArrayOutOfRangeFails)
{
  check_bmc("unpacked_array_out_of_range_fails.sv", 0);
}

// An unpacked array crossing a module boundary: passed into a
// submodule, read there, and driven back out of another. A whole
// array has no bits to splice, so each connection is one term shared
// by both sides.
TEST_P(SVUnitTests, UnpackedArrayPort)
{
  check_prover<KInduction>("unpacked_array_port.sv", 8, ProverResult::TRUE);
}

// Without this, an array the child never actually reached would
// satisfy the holds case above just as well.
TEST_P(SVUnitTests, UnpackedArrayPortFails)
{
  check_bmc("unpacked_array_port_fails.sv", 2);
}

// A negative declared lower bound, a runtime bit position inside a
// runtime-indexed element, and `++`/`--` on an element. The negative
// case asserts cells either side of zero, so an offset computed the
// wrong way round is refuted rather than merely encoded.
TEST_P(SVUnitTests, UnpackedArraySmallForms)
{
  check_prover<KInduction>(
      "unpacked_array_small_forms.sv", 8, ProverResult::TRUE);
}

// More than one dimension: an array of arrays, so each index is a
// Select going down and a Store coming back out. The dimensions are
// unequal and neighbours in both are asserted, so an index applied
// to the wrong dimension is refuted rather than quietly working.
TEST_P(SVUnitTests, UnpackedArrayMultiDim)
{
  // Same pinned-bitwuzla limitation as UnpackedArrayWholeOps: the
  // `'{default: ...}` fills are constant arrays, and comparing one
  // in the same query warns "equality over constant arrays not fully
  // supported yet" and gives up. cvc5 proves it. Drop the special
  // case once bitwuzla is updated.
  ProverResult expected = GetParam() == smt::SolverEnum::BZLA
                              ? ProverResult::UNKNOWN
                              : ProverResult::TRUE;
  check_prover<KInduction>("unpacked_array_2d.sv", 8, expected);
}

// Arrays built outside a clocked block: a combinational lookup
// table, a dynamic-index write over a whole-array default, and a
// whole-array assignment. Writes compose in order and are pinned by
// one constraint per array when the block ends.
TEST_P(SVUnitTests, UnpackedArrayCombWrite)
{
  check_prover<KInduction>("unpacked_array_comb.sv", 6, ProverResult::TRUE);
}

// The same for an `initial` block, which constrains the initial
// state rather than a next one.
TEST_P(SVUnitTests, UnpackedArrayInitialWrite)
{
  check_prover<KInduction>("unpacked_array_initial.sv", 6, ProverResult::TRUE);
}

// Assigning, copying and comparing whole arrays, none of which needs
// the array taken apart -- SMT arrays support all three natively.
TEST_P(SVUnitTests, UnpackedArrayWholeOps)
{
  // The pinned bitwuzla warns "equality over constant arrays not
  // fully supported yet" and gives up when a whole-array comparison
  // meets the constant array the reset assigns; cvc5 proves it. Drop
  // the special case once bitwuzla is updated.
  ProverResult expected = GetParam() == smt::SolverEnum::BZLA
                              ? ProverResult::UNKNOWN
                              : ProverResult::TRUE;
  check_prover<KInduction>("unpacked_array_whole_ops.sv", 12, expected);
}

// ---------------------------------------------------------------------------
// Packed structs
// ---------------------------------------------------------------------------

// Struct-field nonblocking assignment (`p.cnt <= ...`, `s.a.x <= ...`):
// resolve_lvalue() narrows the inner base's bit range by the field's
// own bitOffset, mirroring the read-side MemberAccess case in
// expr_to_term(), so the field is registered as a real assign_next()
// target rather than staying a free state variable.
TEST_P(SVUnitTests, PackedStructFieldState) { check_bmc("struct_state.sv", 5); }

TEST_P(SVUnitTests, PackedStructNested) { check_bmc("struct_nested.sv", 4); }

// Struct-typed *ports* only exercise the read-side MemberAccess path,
// not the field-write path above.
TEST_P(SVUnitTests, TypedefStructPort)
{
  check_bmc("typedef_struct_port.sv", 2);
}

// ---------------------------------------------------------------------------
// Packed enums
// ---------------------------------------------------------------------------

// Referencing one of an enum's own named literals (IDLE/REQ/ACK):
// lookup_symbol() resolves an EnumValueSymbol the same way it resolves
// a ParameterSymbol (both are elaboration-time constants slang has
// already evaluated).
//
// Each FSM fixture below has a `default:` case arm, exercising the
// case-statement default-arm-only-when-no-match behavior tested more
// directly by CaseStatementDefaultOnlyWhenNoMatch in
// test_systemverilog_statements.cpp.
TEST_P(SVUnitTests, PackedEnumStateMachine) { check_bmc("enum_fsm.sv", 3); }

TEST_P(SVUnitTests, PackedEnumStateMachineHolds)
{
  check_bmc("enum_fsm_holds.sv", 5, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, PackedArrayOfEnums) { check_bmc("array_of_enums.sv", 3); }

TEST_P(SVUnitTests, EnumCastFromInt) { check_bmc("enum_cast.sv", 2); }

// ---------------------------------------------------------------------------
// Packed unions
// ---------------------------------------------------------------------------

// Packed-union member access aliases at bit offset 0 for every member
// (unlike struct members, which are packed end-to-end): writing
// through `.b` and reading back through `.parts.hi`/`.parts.lo` is
// bit-consistent.
TEST_P(SVUnitTests, PackedUnionOverlap)
{
  check_bmc("union_overlap.sv", 6, ProverResult::UNKNOWN);
}

// Packed-union construction via `'{default: ...}` (distinct from
// PackedUnionOverlap's member-access above) is ordinary synthesizable
// RTL, not out of scope -- this encoder already supports packed
// unions generally (see PackedUnionOverlap), but
// expr_to_term()'s StructuredAssignmentPattern case only builds a
// PackedStructType target, so a union canonical type throws instead of
// Building a tagged union value, the other half of matching one --
// the same LRM 7.3.2 layout in the other direction, with the bits
// the standard leaves undefined left unconstrained rather than
// zeroed.
TEST_P(SVUnitTests, TaggedUnionConstruction)
{
  check_prover<KInduction>("tagged_union_construct.sv", 6, ProverResult::TRUE);
}

// A `void` member, which exists so the tag can carry everything.
TEST_P(SVUnitTests, TaggedUnionVoidMember)
{
  check_prover<KInduction>("tagged_union_void.sv", 6, ProverResult::TRUE);
}

// The packed qualifier is what makes any of this possible: an
// unpacked union has no required representation, so its tag has no
// position to read.
TEST_P(SVUnitTests, Unsupported_UnpackedTaggedUnion)
{
  expect_encode_throws("unpacked_tagged_union.sv");
}

// enforcing the fixture's own reset-value invariant.
TEST_P(SVUnitTests, Gap_UnionAssignmentPatternLiteral)
{
  check_bmc("union_literal.sv", 4, ProverResult::UNKNOWN);
}

// ---------------------------------------------------------------------------
// typedef: a typedef'd plain vector is just its underlying sort.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, TypedefVectorWidth) { check_bmc("typedef_vector.sv", 5); }

// ---------------------------------------------------------------------------
// 2-state (`bit`) vs 4-state (`logic`) parity: Pono's SMT bitvector
// model has no X/Z state, so `bit` and `logic` registers updated
// identically are numerically indistinguishable (checked as both a
// holding and a violated property below).
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, BitVsLogicParityHolds)
{
  check_bmc("bit_vs_logic_parity.sv", 10, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, BitVsLogicParityMismatchFails)
{
  check_bmc("bit_vs_logic_parity_mismatch.sv", 1);
}

// ---------------------------------------------------------------------------
// `signed` types: width-extension, comparison, and division must all use
// two's-complement signed semantics, not treat the raw bit pattern as
// unsigned.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, SignedSignExtend)
{
  check_bmc("signed_sign_extend.sv", 0, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, SignedCompare)
{
  check_bmc("signed_compare.sv", 0, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, SignedDivide)
{
  check_bmc("signed_divide.sv", 0, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, SignedCast)
{
  check_bmc("signed_cast.sv", 0, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, SignedContextWidth)
{
  check_bmc("signed_context_width.sv", 0, ProverResult::UNKNOWN);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVTypesTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

// A constant pattern lists its values in declared-index order, which
// runs the opposite way from the normalized index a descending range
// uses -- so filling by position reverses the array.

TEST_P(SVUnitTests, UnpackedArrayDescendingPattern)
{
  check_prover<KInduction>(
      "unpacked_array_descending.sv", 10, ProverResult::TRUE);
}

TEST_P(SVUnitTests, UnpackedArrayDescendingPatternFails)
{
  check_bmc("unpacked_array_descending_fails.sv", 2);
}

// A pattern whose values are not elaboration-time constants: one
// term per element, built as stores over a seed rather than folded
// into a single constant array.

TEST_P(SVUnitTests, UnpackedPatternNonConstant)
{
  check_prover<KInduction>(
      "unpacked_pattern_nonconstant.sv", 10, ProverResult::TRUE);
}

TEST_P(SVUnitTests, UnpackedPatternNonConstantFails)
{
  check_bmc("unpacked_pattern_nonconstant_fails.sv", 2);
}

// Unpacked-array nets. Driving one is a constraint rather than an
// assignment, and what is left undriven stays free -- which is what
// an undriven net is.

TEST_P(SVUnitTests, UnpackedArrayNet)
{
  check_prover<KInduction>("unpacked_array_net.sv", 8, ProverResult::TRUE);
}

TEST_P(SVUnitTests, UnpackedArrayNetUndrivenIsFree)
{
  check_bmc("unpacked_array_net_fails.sv", 1);
}

TEST_P(SVUnitTests, UnpackedArrayOutPort)
{
  check_prover<KInduction>("unpacked_array_out_port.sv", 8, ProverResult::TRUE);
}

// Unpacked structs. No bit width of their own, but a selectable
// one -- the space each field's bitOffset is measured in -- so a
// flat layout there makes a field a bit range, as for a packed
// struct. Covers a plain variable, an array element, and a whole
// copy; the fields differ in width so an overlapping offset shows.

TEST_P(SVUnitTests, UnpackedStruct)
{
  check_prover<KInduction>("unpacked_struct.sv", 10, ProverResult::TRUE);
}

TEST_P(SVUnitTests, UnpackedStructFieldsDoNotOverlap)
{
  check_bmc("unpacked_struct_fails.sv", 2);
}

// An unpacked array connected to a slice of a parent array: the two
// have different lengths, so the port gets its own array tied to
// the parent's element by element. (A *concatenation* target is not
// legal SystemVerilog for an unpacked array at all -- slang rejects
// it as not assignable.)

TEST_P(SVUnitTests, UnpackedArrayPortSlice)
{
  check_prover<KInduction>(
      "unpacked_array_port_slice.sv", 8, ProverResult::TRUE);
}

TEST_P(SVUnitTests, UnpackedArrayPortSliceOffsetIsReal)
{
  check_bmc("unpacked_array_port_slice_fails.sv", 1);
}

}  // namespace pono_tests
