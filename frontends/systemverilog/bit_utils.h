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
class Type;
}  // namespace slang::ast

namespace pono {

/** Convert an integral slang type to a BV sort (even for a 1-bit type).
 *  @param solver the solver to build the sort with
 *  @param type the slang type
 *  @return the corresponding BV sort; throws for a non-integral type
 */
smt::Sort type_to_sort(const smt::SmtSolver & solver,
                       const slang::ast::Type & type);

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

/** Build a partial-write term like replace_bits(), but for a
 *  runtime-variable element index (`arr[idx] = slice` where `idx` is not
 *  a compile-time constant): shifts an `elem_w`-wide window of ones into
 *  position `idx * elem_w` and uses it to mask the shifted-into-position
 *  `slice` into `base`, leaving every other element unchanged. Mirrors
 *  the shift+extract technique used for dynamic-index *reads*.
 */
smt::Term replace_bits_dynamic(const smt::SmtSolver & solver,
                               const smt::Term & base,
                               const smt::Term & slice,
                               const smt::Term & idx,
                               uint64_t elem_w);

}  // namespace pono
