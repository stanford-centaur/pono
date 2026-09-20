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
  if (!arr.elementType.isIntegral()) {
    // A packed struct or union element is integral and so does pass:
    // what this rules out is a further unpacked dimension, or an
    // unpacked struct/union element.
    throw PonoException(
        "SystemVerilogEncoder: only unpacked arrays of an integral element "
        "type are supported, so not '"
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

  throw PonoException("SystemVerilogEncoder: unsupported type kind");
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

Term replace_bits_dynamic(const SmtSolver & solver,
                          const Term & base,
                          const Term & slice,
                          const Term & idx,
                          uint64_t elem_w,
                          uint64_t base_offset)
{
  // Index arithmetic and bit masks below -- zero-extend throughout;
  // `slice` is padded before shifting into position, but `mask`
  // clears every bit outside the shifted elem_w-wide window
  // regardless, so the padding bits of `slice` never affect the
  // result either way.
  require_bv(base, "replace_bits_dynamic()");
  uint64_t base_w = base->get_sort()->get_width();
  Sort base_sort = solver->make_sort(BV, base_w);
  Term idx_ext = resize_to(solver, idx, base_w, false);
  Term shift_amount = idx_ext;
  if (elem_w != 1) {
    Term elem_w_term = solver->make_term(elem_w, base_sort);
    shift_amount = solver->make_term(BVMul, idx_ext, elem_w_term);
  }
  if (base_offset != 0) {
    shift_amount = solver->make_term(
        BVAdd, shift_amount, solver->make_term(base_offset, base_sort));
  }
  Term elem_ones = solver->make_term(
      BVNot, solver->make_term(0, solver->make_sort(BV, elem_w)));
  Term mask = solver->make_term(
      BVShl, resize_to(solver, elem_ones, base_w, false), shift_amount);
  Term slice_shifted = solver->make_term(
      BVShl, resize_to(solver, slice, base_w, false), shift_amount);
  Term cleared = solver->make_term(BVAnd, base, solver->make_term(BVNot, mask));
  return solver->make_term(
      BVOr, cleared, solver->make_term(BVAnd, slice_shifted, mask));
}

}  // namespace pono
