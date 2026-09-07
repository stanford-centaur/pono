/*!
 * \file bit_utils.cpp
 * \brief Pure bit/type-manipulation helpers shared across the SV encoder.
 * \author Áron Ricardo Perez-Lopez
 * \date 2026
 * \copyright See the LICENSE file in the top-level source directory.
 */
#include "frontends/systemverilog/bit_utils.h"

#include <vector>

#include "slang/ast/types/Type.h"
#include "utils/exceptions.h"

using namespace smt;
using namespace std;

namespace pono {

Sort type_to_sort(const SmtSolver & solver, const slang::ast::Type & type)
{
  if (type.isIntegral()) {
    uint64_t width = type.getBitWidth();
    if (width == 0) {
      throw PonoException("SystemVerilogEncoder: zero-width integral type");
    }
    return solver->make_sort(BV, width);
  }

  throw PonoException("SystemVerilogEncoder: unsupported type kind");
}

Term slice_bits(const SmtSolver & solver,
                const Term & base,
                uint64_t lo,
                uint64_t hi)
{
  if (!base) return Term();
  uint64_t w = base->get_sort()->get_width();
  if (lo == 0 && hi == w - 1) return base;
  return solver->make_term(Op(Extract, hi, lo), base);
}

Term resize_to(const SmtSolver & solver,
               const Term & t,
               uint64_t target_width,
               bool is_signed)
{
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
                          uint64_t elem_w)
{
  // Index arithmetic and bit masks below -- zero-extend throughout;
  // `slice` is padded before shifting into position, but `mask`
  // clears every bit outside the shifted elem_w-wide window
  // regardless, so the padding bits of `slice` never affect the
  // result either way.
  uint64_t base_w = base->get_sort()->get_width();
  Sort base_sort = solver->make_sort(BV, base_w);
  Term idx_ext = resize_to(solver, idx, base_w, false);
  Term shift_amount = idx_ext;
  if (elem_w != 1) {
    Term elem_w_term = solver->make_term(elem_w, base_sort);
    shift_amount = solver->make_term(BVMul, idx_ext, elem_w_term);
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
