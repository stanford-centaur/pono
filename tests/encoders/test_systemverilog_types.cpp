#ifdef WITH_SLANG

#include "sv_test_fixture.h"

using namespace pono;
using namespace smt;

namespace pono_tests {

// ---------------------------------------------------------------------------
// Packed multi-dimensional arrays (array-of-vectors) -- these pass: both
// constant- and variable-index reads/writes on a packed array are already
// supported.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, PackedArrayConstIndexChain)
{
  check_bmc("array_const_index.sv", 7);
}

TEST_P(SVUnitTests, PackedArrayDynIndexReadWrite)
{
  check_bmc("array_dyn_index.sv", 16);
}

// ---------------------------------------------------------------------------
// Packed structs
// ---------------------------------------------------------------------------

// FIXED: find_lhs_base()/resolve_lvalue() in systemverilog_encoder.cpp
// used to have no ExpressionKind::MemberAccess case, so a struct-field
// nonblocking assignment (`p.cnt <= ...`, `s.a.x <= ...`) failed to
// resolve as an lvalue and was silently dropped -- the field was never
// registered as an assign_next() target and stayed a completely free
// state variable every cycle. Both helpers now narrow the inner base's
// bit range by the field's own bitOffset, the same way the read-side
// MemberAccess case in expr_to_term() already did.
TEST_P(SVUnitTests, PackedStructFieldState) { check_bmc("struct_state.sv", 5); }

TEST_P(SVUnitTests, PackedStructNested) { check_bmc("struct_nested.sv", 4); }

// Passes: struct-typed *ports* only exercise the (supported) read-side
// MemberAccess path, not the field-write path fixed above.
TEST_P(SVUnitTests, TypedefStructPort)
{
  check_bmc("typedef_struct_port.sv", 2);
}

// ---------------------------------------------------------------------------
// Packed enums
// ---------------------------------------------------------------------------

// FIXED: lookup_symbol() had no path for an EnumValueSymbol -- only
// declare_variables()-tracked registers/wires/parameters were ever
// inserted into symbol_to_term_ -- so referencing one of an enum's own
// named literals (IDLE/REQ/ACK, used exactly as ordinary SV permits
// once the enum is declared) threw "unknown symbol".  Declaring an
// enum-typed *variable* already worked fine (its type is integral, so
// type_to_sort() succeeds); lookup_symbol() now also resolves an
// EnumValueSymbol the same way it already resolved a ParameterSymbol
// (both are elaboration-time constants slang has already evaluated).
//
// Each FSM fixture below has a `default:` case arm -- exercising one
// was needed to fully test these fixtures, and along the way it
// surfaced (and, separately, fixed -- see
// CaseStatementDefaultOnlyWhenNoMatch in test_systemverilog_statements
// .cpp) a distinct, pre-existing bug where a `case` statement's
// `default:` arm applied unconditionally alongside whichever other
// item already matched, sticking every such FSM at its default value
// forever.
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

// Passes: packed-union member access correctly aliases at bit offset
// 0 for every member (unlike struct members, which are packed
// end-to-end) -- writing through `.b` and reading back through
// `.parts.hi`/`.parts.lo` is bit-consistent.
TEST_P(SVUnitTests, PackedUnionOverlap)
{
  check_bmc("union_overlap.sv", 6, ProverResult::UNKNOWN);
}

// ---------------------------------------------------------------------------
// typedef -- passes: a typedef'd plain vector is just its underlying sort.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, TypedefVectorWidth) { check_bmc("typedef_vector.sv", 5); }

// ---------------------------------------------------------------------------
// 2-state (`bit`) vs 4-state (`logic`) parity -- passes both directions:
// Pono's SMT bitvector model has no X/Z state, so `bit` and `logic`
// registers updated identically are numerically indistinguishable.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, BitVsLogicParityHolds)
{
  check_bmc("bit_vs_logic_parity.sv", 10, ProverResult::UNKNOWN);
}

TEST_P(SVUnitTests, BitVsLogicParityMismatchFails)
{
  check_bmc("bit_vs_logic_parity_mismatch.sv", 1);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVTypesTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests

#endif
