/*!
 * \file bit_utils.h
 * \brief Pure bit/type-manipulation helpers shared across the SV encoder.
 * \author Áron Ricardo Perez-Lopez
 * \date 2026
 * \copyright See the LICENSE file in the top-level source directory.
 *
 * Free functions with no state of their own -- each takes the solver
 * explicitly and returns a new Term/Sort, so any class in this directory
 * can call them without depending on anything but smt-switch.
 */
#pragma once

#include <cstdint>

#include "smt-switch/smt.h"

namespace slang::ast {
class FixedSizeUnpackedArrayType;
class PackedArrayType;
class Type;
}  // namespace slang::ast

namespace pono {

/** How a fixed-size unpacked array's declared SV index range maps onto
 *  the index sort of the SMT array modelling it.
 *
 *  Element 0 of the SMT array is the element declared at `lower`, so
 *  `mem[3:18]` reads `mem[i]` as `Select(mem, i - 3)`.  Every
 *  Select/Store has to apply the same mapping -- getting it wrong in
 *  one place is silently wrong rather than broken -- so go through
 *  normalize_array_index() rather than re-deriving it.
 */
struct UnpackedArrayInfo
{
  uint64_t depth;          ///< number of elements
  int64_t lower;           ///< SV index of element 0, which may be negative
  uint64_t index_width;    ///< width of the SMT array's index sort
  smt::Sort element_sort;  ///< sort of one element
};

/** Describe `arr`'s index mapping.  Callers narrow the type first --
 *  `type.getCanonicalType().kind == SymbolKind::FixedSizeUnpackedArrayType`
 *  then `.as<FixedSizeUnpackedArrayType>()` -- so there is nothing here
 *  to report absence about.
 *
 *  Throws only for an element type with no sort of its own, such as
 *  an unpacked struct. A negative declared lower bound, a further
 *  dimension, and a packed struct or union element are all modelled.
 */
UnpackedArrayInfo unpacked_array_info(
    const smt::SmtSolver & solver,
    const slang::ast::FixedSizeUnpackedArrayType & arr);

/** Map an SV element index onto the array's SMT index sort: subtract
 *  the declared lower bound, then resize to the index width.
 *
 *  Also reports, through `in_range` when non-null, whether the index
 *  names a cell the array actually has. Left untouched when no index
 *  can miss -- a zero-based array whose depth fills its address
 *  space, read by an index no wider than that space. Out-of-range
 *  accesses must consult it: the truncation to `index_width` would
 *  otherwise silently fold, say, `m[23]` of a `m[3:18]` onto `m[7]`. */
smt::Term normalize_array_index(const smt::SmtSolver & solver,
                                const smt::Term & idx,
                                const UnpackedArrayInfo & info,
                                smt::Term * in_range = nullptr);

/** Convert a slang type to an SMT sort: a BV for an integral type (even
 *  a 1-bit one), an ARRAY for a fixed-size unpacked array of an
 *  integral element type.
 *  @param solver the solver to build the sort with
 *  @param type the slang type
 *  @return the corresponding sort; throws for anything else
 */
smt::Sort type_to_sort(const smt::SmtSolver & solver,
                       const slang::ast::Type & type);

/** Throw unless `t` is bit-vector-sorted, naming `who` and the sort.
 *
 *  `get_sort()->get_width()` is the backend's bit-vector-size query,
 *  so reaching it with an array- or Bool-sorted term aborts inside
 *  cvc5/bitwuzla rather than raising a PonoException.  Anything about
 *  to do width arithmetic on a term that could be an unpacked array --
 *  the bit helpers below, and the operand preparation in
 *  expr_encoder.cpp -- checks first.
 */
void require_bv(const smt::Term & t, const char * who);

/** Turn a packed array's declared SV index into an ordinal counted
 *  from its least significant element, which is what a bit offset is
 *  measured in.
 *
 *  Only `[n:0]` makes the two the same. A range that starts
 *  elsewhere is shifted by its lower bound, and an ascending one
 *  (`[0:n]`, whose leftmost element is the most significant) runs the
 *  other way entirely -- so using the index itself reads and writes
 *  the wrong element rather than failing.
 *
 *  @return false if `idx` is outside the declared range, leaving
 *          `ordinal` untouched; the caller decides whether that is an
 *          X-valued read or an error.
 */
bool packed_element_ordinal(const slang::ast::PackedArrayType & arr,
                            int64_t idx,
                            uint64_t & ordinal);

/** How many bits a value of this type occupies, which is not always
 *  its bit width: an unpacked struct has none of its own, but it has
 *  a *selectable* width -- the space its fields' own bitOffsets are
 *  measured in, and so the one a field's bit range lives in.
 *
 *  Use this wherever the question is how wide a variable, field,
 *  element or port is. `getBitWidth()` still answers for an
 *  expression's result, a literal or a loop counter, all of which
 *  are integral by construction.
 */
uint64_t value_width(const slang::ast::Type & type);

/** Re-order a bit stream into `slice`-wide blocks, as a `<<`
 *  streaming concatenation does (LRM 11.4.14.2). Blocks are cut from
 *  the right-hand end, so when `slice` does not divide the width it
 *  is the *left-most* block that is short; the block order is then
 *  reversed, leaving that short block at the right-hand end.
 *
 *  `slice` is zero for `>>`, which re-orders nothing, and so is a
 *  no-op here.
 */
smt::Term stream_reorder(const smt::SmtSolver & solver,
                         const smt::Term & value,
                         uint64_t slice);

/** Undo stream_reorder(), recovering the generic stream that a `<<`
 *  produced -- the operation a streaming concatenation performs when
 *  it is an assignment target (LRM 11.4.14.3).
 *
 *  The two are the same permutation only when `slice` divides the
 *  width. Otherwise the short block sits at the opposite end, so
 *  this cuts it from the right-hand end first and moves it to the
 *  left.
 */
smt::Term stream_unreorder(const smt::SmtSolver & solver,
                           const smt::Term & value,
                           uint64_t slice);

/** Extract bits [lo, hi] from `base`, or return `base` unchanged when
 *  [lo, hi] already covers its whole width. Returns a null Term if
 *  `base` is null.
 */
smt::Term slice_bits(const smt::SmtSolver & solver,
                     const smt::Term & base,
                     uint64_t lo,
                     uint64_t hi);

/** Ensure a term has the expected bit-width, extending or truncating as
 *  needed.
 *  @param solver the solver to build the resized term with
 *  @param t the term to resize
 *  @param target_width the desired width
 *  @param is_signed if the term grows, sign-extend (replicate the top
 *         bit) instead of zero-extending -- use this whenever `t` is the
 *         value of a `signed`-typed SystemVerilog expression being
 *         widened; every other use (bit-range bookkeeping, unsigned
 *         types, etc.) should pass false. Required (no default) so every
 *         call site makes an explicit, reviewed choice instead of
 *         silently defaulting to zero-extension.
 *  @return the resized term
 */
smt::Term resize_to(const smt::SmtSolver & solver,
                    const smt::Term & t,
                    uint64_t target_width,
                    bool is_signed);

/** Build a partial-write term: take the full-width `base` and return a
 *  term equal to `base` everywhere except bits [lo .. hi], which take
 *  their values from `slice`. A full-width write (lo == 0 && hi ==
 *  width(base)-1) is just `slice`.
 */
smt::Term replace_bits(const smt::SmtSolver & solver,
                       const smt::Term & base,
                       const smt::Term & slice,
                       uint64_t lo,
                       uint64_t hi);

/** Splice `width` bits of `slice` into `base` at a runtime bit
 *  `position`, leaving every other bit alone. The element-indexed
 *  form below is this with the position worked out from an index;
 *  callers that already know the position in bits -- because it was
 *  rebased onto some other signal -- want this one.
 */
smt::Term replace_bits_at(const smt::SmtSolver & solver,
                          const smt::Term & base,
                          const smt::Term & slice,
                          const smt::Term & position,
                          uint64_t width);

/** Build a partial-write term like replace_bits(), but for a
 *  runtime-variable element index (`arr[idx] = slice` where `idx` is not
 *  a compile-time constant): shifts an `elem_w`-wide window of ones into
 *  position `base_offset + idx * elem_w` and uses it to mask the
 *  shifted-into-position `slice` into `base`, leaving every other
 *  element unchanged. Mirrors the shift+extract technique used for
 *  dynamic-index *reads*.
 *
 *  `base_offset` is where the indexed range itself starts within
 *  `base`, which is nonzero when the select sits on a sub-range rather
 *  than on the whole variable (`p[2][idx]`, `s.field[idx]`).
 */
smt::Term replace_bits_dynamic(const smt::SmtSolver & solver,
                               const smt::Term & base,
                               const smt::Term & slice,
                               const smt::Term & idx,
                               uint64_t elem_w,
                               uint64_t base_offset = 0);

}  // namespace pono
