/*!
 * \file bit_utils.cpp
 * \brief Pure bit/type-manipulation helpers shared across the SV encoder.
 * \author Áron Ricardo Perez-Lopez
 * \date 2026
 * \copyright See the LICENSE file in the top-level source directory.
 */
#include "frontends/systemverilog/bit_utils.h"

#include <algorithm>
#include <vector>

#include "slang/ast/symbols/VariableSymbols.h"
#include "slang/ast/types/AllTypes.h"
#include "slang/ast/types/Type.h"
#include "utils/exceptions.h"

using namespace smt;
using namespace std;

namespace pono {

namespace {

/** Narrowest index width that can address `depth` elements. */
uint64_t index_width_for_depth(uint64_t depth)
{
  uint64_t w = 1;
  while ((uint64_t{ 1 } << w) < depth) ++w;
  return w;
}

}  // namespace

void require_bv(const Term & t, const char * who)
{
  if (t->get_sort()->get_sort_kind() != BV) {
    throw PonoException("SystemVerilogEncoder: " + std::string(who)
                        + " expects a bit-vector, got sort "
                        + t->get_sort()->to_string());
  }
}

UnpackedArrayInfo unpacked_array_info(
    const SmtSolver & solver,
    const slang::ast::FixedSizeUnpackedArrayType & arr)
{
  // A further unpacked dimension is fine -- type_to_sort() recurses
  // and builds an array of arrays. A packed struct or union element
  // is integral and passes too, and an unpacked struct has a flat
  // layout of its own. What this rules out is an element with no
  // sort at all; type_to_sort() reports which, since it is the one
  // that knows why.
  auto elem_kind = arr.elementType.getCanonicalType().kind;
  bool composite =
      elem_kind == slang::ast::SymbolKind::FixedSizeUnpackedArrayType
      || elem_kind == slang::ast::SymbolKind::UnpackedStructType;
  if (!arr.elementType.isIntegral() && !composite) {
    throw PonoException(
        "SystemVerilogEncoder: an unpacked array's element type must be "
        "integral, an unpacked struct, or another unpacked array, so not '"
        + std::string(arr.elementType.toString()) + "'");
  }

  uint64_t depth = arr.range.fullWidth();
  return UnpackedArrayInfo{ depth,
                            static_cast<int64_t>(arr.range.lower()),
                            index_width_for_depth(depth),
                            type_to_sort(solver, arr.elementType) };
}

Term normalize_array_index(const SmtSolver & solver,
                           const Term & idx,
                           const UnpackedArrayInfo & info,
                           Term * in_range)
{
  require_bv(idx, "normalize_array_index()");
  uint64_t idx_w = idx->get_sort()->get_width();

  // Can the source index name something outside the declared range?
  // Only then is a check worth building: a zero-based array whose
  // depth fills its address space cannot be missed by an index no
  // wider than that space.
  bool checkable = info.lower != 0 || idx_w > info.index_width
                   || info.depth != (uint64_t{ 1 } << info.index_width);

  if (info.lower == 0 && !checkable) {
    return resize_to(solver, idx, info.index_width, /*is_signed=*/false);
  }

  // Work wide enough that the subtraction below cannot wrap for any
  // index the array actually has, so an index outside the range
  // underflows to a value the check still rejects and a valid one
  // survives the truncation.
  uint64_t magnitude = info.lower < 0 ? static_cast<uint64_t>(-info.lower)
                                      : static_cast<uint64_t>(info.lower);
  uint64_t w = std::max<uint64_t>(idx_w, info.index_width) + 1;
  while ((uint64_t{ 1 } << (w - 1)) < magnitude + info.depth) ++w;
  Sort wide_sort = solver->make_sort(BV, w);
  // A negative lower bound means the declared range includes negative
  // indices, so the index expression is a signed value and has to
  // widen as one.
  Term wide = resize_to(solver, idx, w, /*is_signed=*/info.lower < 0);
  Term lower_term = solver->make_term(magnitude, wide_sort);
  if (info.lower < 0) lower_term = solver->make_term(BVNeg, lower_term);
  Term offset =
      info.lower == 0 ? wide : solver->make_term(BVSub, wide, lower_term);
  if (in_range && checkable) {
    *in_range = solver->make_term(
        BVUlt, offset, solver->make_term(info.depth, wide_sort));
  }
  return resize_to(solver, offset, info.index_width, /*is_signed=*/false);
}

Sort type_to_sort(const SmtSolver & solver, const slang::ast::Type & type)
{
  if (type.isIntegral()) {
    uint64_t width = type.getBitWidth();
    if (width == 0) {
      throw PonoException("SystemVerilogEncoder: zero-width integral type");
    }
    return solver->make_sort(BV, width);
  }

  const slang::ast::Type & ct = type.getCanonicalType();
  if (ct.kind == slang::ast::SymbolKind::FixedSizeUnpackedArrayType) {
    UnpackedArrayInfo info = unpacked_array_info(
        solver, ct.as<slang::ast::FixedSizeUnpackedArrayType>());
    return solver->make_sort(
        ARRAY, solver->make_sort(BV, info.index_width), info.element_sort);
  }

  if (ct.kind == slang::ast::SymbolKind::UnpackedStructType) {
    // An unpacked struct has no bit width of its own, but it does
    // have a selectable one -- the space its fields' own bitOffsets
    // are measured in. Laying it out flat in that space is what
    // makes a field a bit range, which is how a packed struct is
    // already read and written.
    auto & st = ct.as<slang::ast::UnpackedStructType>();
    for (auto * field : st.fields) {
      if (!field->getType().isIntegral()) {
        throw PonoException(
            "SystemVerilogEncoder: field '" + std::string(field->name)
            + "' of an unpacked struct is not integral, so the struct has "
              "no flat layout to give it a bit range in");
      }
    }
    if (st.selectableWidth == 0) {
      throw PonoException("SystemVerilogEncoder: zero-width unpacked struct");
    }
    return solver->make_sort(BV, st.selectableWidth);
  }

  throw PonoException("SystemVerilogEncoder: unsupported type kind");
}

uint64_t value_width(const slang::ast::Type & type)
{
  const slang::ast::Type & ct = type.getCanonicalType();
  if (ct.kind == slang::ast::SymbolKind::UnpackedStructType) {
    return ct.as<slang::ast::UnpackedStructType>().selectableWidth;
  }
  return type.getBitWidth();
}

bool packed_element_ordinal(const slang::ast::PackedArrayType & arr,
                            int64_t idx,
                            uint64_t & ordinal)
{
  int64_t lower = arr.range.lower();
  int64_t upper = arr.range.upper();
  if (idx < lower || idx > upper) return false;
  // The left bound is the most significant end either way, so a
  // descending range counts up from its lower bound and an ascending
  // one counts down from its upper.
  bool descending = arr.range.left >= arr.range.right;
  ordinal = static_cast<uint64_t>(descending ? idx - lower : upper - idx);
  return true;
}

namespace {

// The block partition both re-orderings work on: `q` whole blocks and
// a remainder of `rem` bits. Reports false when there is nothing to
// re-order -- a `>>` stream, or a single block covering everything.
bool stream_blocks(const Term & value,
                   uint64_t slice,
                   uint64_t & width,
                   uint64_t & q,
                   uint64_t & rem)
{
  require_bv(value, "stream re-ordering");
  width = value->get_sort()->get_width();
  if (slice == 0 || slice >= width) return false;
  rem = width % slice;
  q = (width - rem) / slice;
  return q > 1 || rem != 0;
}

Term concat_all(const SmtSolver & solver, const std::vector<Term> & pieces)
{
  Term result;
  for (const Term & p : pieces) {
    result = result ? solver->make_term(Concat, result, p) : p;
  }
  return result;
}

}  // namespace

TaggedUnionLayout tagged_union_layout(const slang::ast::Type & type,
                                      const char * what)
{
  using namespace slang::ast;
  const Type & canon = type.getCanonicalType();
  if (canon.kind != SymbolKind::PackedUnionType
      || !canon.as<PackedUnionType>().isTagged) {
    throw PonoException(std::string("SystemVerilogEncoder: ") + what
                        + " needs a packed tagged union, and '"
                        + type.toString()
                        + "' is not one -- an unpacked union has no required "
                          "representation, so its tag has no position to read");
  }
  auto & pu = canon.as<PackedUnionType>();
  return { pu.getBitWidth(), pu.tagBits };
}

Term stream_reorder(const SmtSolver & solver,
                    const Term & value,
                    uint64_t slice)
{
  uint64_t width = 0, q = 0, rem = 0;
  if (!stream_blocks(value, slice, width, q, rem)) return value;
  // Blocks numbered from the right of the input come out left-first,
  // and the short left-over -- the input's top `rem` bits -- last.
  std::vector<Term> pieces;
  pieces.reserve(q + 1);
  for (uint64_t i = 0; i < q; i++) {
    pieces.push_back(
        solver->make_term(Op(Extract, (i + 1) * slice - 1, i * slice), value));
  }
  if (rem) {
    pieces.push_back(
        solver->make_term(Op(Extract, width - 1, q * slice), value));
  }
  return concat_all(solver, pieces);
}

Term stream_unreorder(const SmtSolver & solver,
                      const Term & value,
                      uint64_t slice)
{
  uint64_t width = 0, q = 0, rem = 0;
  if (!stream_blocks(value, slice, width, q, rem)) return value;
  // The mirror image: the short block is the input's *bottom* `rem`
  // bits and belongs at the top, ahead of the whole blocks taken
  // right to left.
  std::vector<Term> pieces;
  pieces.reserve(q + 1);
  if (rem) {
    pieces.push_back(solver->make_term(Op(Extract, rem - 1, 0), value));
  }
  for (uint64_t i = 0; i < q; i++) {
    pieces.push_back(solver->make_term(
        Op(Extract, rem + (i + 1) * slice - 1, rem + i * slice), value));
  }
  return concat_all(solver, pieces);
}

Term slice_bits(const SmtSolver & solver,
                const Term & base,
                uint64_t lo,
                uint64_t hi)
{
  if (!base) return Term();
  require_bv(base, "slice_bits()");
  uint64_t w = base->get_sort()->get_width();
  if (lo == 0 && hi == w - 1) return base;
  return solver->make_term(Op(Extract, hi, lo), base);
}

Term resize_to(const SmtSolver & solver,
               const Term & t,
               uint64_t target_width,
               bool is_signed)
{
  require_bv(t, "resize_to()");
  uint64_t current_width = t->get_sort()->get_width();
  if (current_width == target_width) {
    return t;
  }
  if (current_width < target_width) {
    Op ext_op(is_signed ? Sign_Extend : Zero_Extend,
              target_width - current_width);
    return solver->make_term(ext_op, t);
  }
  // Truncate (extract lower bits).
  return solver->make_term(Op(Extract, target_width - 1, 0), t);
}

Term replace_bits(const SmtSolver & solver,
                  const Term & base,
                  const Term & slice,
                  uint64_t lo,
                  uint64_t hi)
{
  require_bv(base, "replace_bits()");
  uint64_t base_w = base->get_sort()->get_width();
  if (lo == 0 && hi == base_w - 1) return slice;
  std::vector<Term> parts;
  if (hi + 1 < base_w) {
    parts.push_back(solver->make_term(Op(Extract, base_w - 1, hi + 1), base));
  }
  parts.push_back(slice);
  if (lo > 0) {
    parts.push_back(solver->make_term(Op(Extract, lo - 1, 0), base));
  }
  Term result = parts[0];
  for (size_t i = 1; i < parts.size(); ++i) {
    result = solver->make_term(Concat, result, parts[i]);
  }
  return result;
}

Term replace_bits_at(const SmtSolver & solver,
                     const Term & base,
                     const Term & slice,
                     const Term & position,
                     uint64_t width)
{
  // Bit masks throughout -- zero-extend everything; `slice` is padded
  // before shifting into position, but `mask` clears every bit
  // outside the shifted width-wide window regardless, so the padding
  // bits never affect the result either way.
  require_bv(base, "replace_bits_at()");
  uint64_t base_w = base->get_sort()->get_width();
  Term shift = resize_to(solver, position, base_w, false);
  Term ones = solver->make_term(
      BVNot, solver->make_term(0, solver->make_sort(BV, width)));
  Term mask =
      solver->make_term(BVShl, resize_to(solver, ones, base_w, false), shift);
  Term shifted =
      solver->make_term(BVShl, resize_to(solver, slice, base_w, false), shift);
  Term cleared = solver->make_term(BVAnd, base, solver->make_term(BVNot, mask));
  return solver->make_term(
      BVOr, cleared, solver->make_term(BVAnd, shifted, mask));
}

Term replace_bits_dynamic(const SmtSolver & solver,
                          const Term & base,
                          const Term & slice,
                          const Term & idx,
                          uint64_t elem_w,
                          uint64_t base_offset)
{
  require_bv(base, "replace_bits_dynamic()");
  uint64_t base_w = base->get_sort()->get_width();
  Sort base_sort = solver->make_sort(BV, base_w);
  Term shift = resize_to(solver, idx, base_w, false);
  if (elem_w != 1) {
    shift =
        solver->make_term(BVMul, shift, solver->make_term(elem_w, base_sort));
  }
  if (base_offset != 0) {
    shift = solver->make_term(
        BVAdd, shift, solver->make_term(base_offset, base_sort));
  }
  return replace_bits_at(solver, base, slice, shift, elem_w);
}

}  // namespace pono
