/*! \file terms.cpp
 *  \brief SystemVerilogEncoder's low-level bit/term helper substrate:
 *         type-to-sort conversion, output-port-alias resolution, symbol
 *         lookup/naming, and the resize/replace-bits term-manipulation
 *         primitives used throughout the rest of the encoder.
 */
#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

#include "frontends/systemverilog/encoder.h"
#include "slang/ast/Symbol.h"
#include "slang/ast/symbols/ParameterSymbols.h"
#include "slang/ast/symbols/ValueSymbol.h"
#include "slang/ast/types/AllTypes.h"
#include "slang/ast/types/Type.h"
#include "smt-switch/smt.h"
#include "utils/exceptions.h"

using namespace smt;
using namespace std;

namespace pono {

Sort SystemVerilogEncoder::type_to_sort(const slang::ast::Type & type)
{
  if (type.isIntegral()) {
    uint64_t width = type.getBitWidth();
    if (width == 0) {
      throw PonoException("SystemVerilogEncoder: zero-width integral type");
    }
    return solver_->make_sort(BV, width);
  }

  throw PonoException("SystemVerilogEncoder: unsupported type kind");
}

// ============================================================================
// Helpers
// ============================================================================

std::vector<SystemVerilogEncoder::ResolvedAliasPiece>
SystemVerilogEncoder::resolve_output_alias_pieces(
    const slang::ast::Symbol * sym,
    uint64_t lo,
    uint64_t hi,
    uint64_t rhs_base) const
{
  auto alias_it = port_output_aliases_.find(sym);
  if (alias_it == port_output_aliases_.end()) {
    return { { sym, lo, hi, rhs_base, rhs_base + (hi - lo) } };
  }
  std::vector<ResolvedAliasPiece> result;
  for (auto & seg : alias_it->second) {
    // Intersect the caller's [lo, hi] window (in sym's own numbering)
    // with this segment's own [port_lo, port_hi] coverage.
    uint64_t ilo = std::max(lo, seg.port_lo);
    uint64_t ihi = std::min(hi, seg.port_hi);
    if (ilo > ihi) continue;
    uint64_t offset = ilo - seg.port_lo;
    uint64_t span = ihi - ilo;
    uint64_t tlo = seg.target_lo + offset;
    uint64_t thi = tlo + span;
    // Recurse in case `seg.target` is itself an output-port alias
    // (e.g. a nested/chained instantiation); rhs_base advances by
    // however far into the caller's own window this segment starts.
    auto sub = resolve_output_alias_pieces(
        seg.target, tlo, thi, rhs_base + (ilo - lo));
    result.insert(result.end(), sub.begin(), sub.end());
  }
  return result;
}

smt::Term SystemVerilogEncoder::slice_bits(const smt::Term & base,
                                           uint64_t lo,
                                           uint64_t hi) const
{
  if (!base) return Term();
  uint64_t w = base->get_sort()->get_width();
  if (lo == 0 && hi == w - 1) return base;
  return solver_->make_term(Op(Extract, hi, lo), base);
}

bool SystemVerilogEncoder::resolve_wire_on_demand(
    const slang::ast::Symbol * sym)
{
  auto it = wire_drivers_.find(sym);
  if (it == wire_drivers_.end()) return false;
  const WireDriver & drv = it->second;
  const void * stmt_key = drv.ca ? static_cast<const void *>(drv.ca)
                                 : static_cast<const void *>(drv.comb);
  if (processed_drivers_.count(stmt_key)) return true;

  if (!resolving_wires_.insert(sym).second) {
    throw PonoException(
        "SystemVerilogEncoder: combinational loop detected involving '"
        + std::string(sym->name) + "'");
  }
  string saved_prefix = prefix_;
  string saved_parent_prefix = parent_prefix_;
  prefix_ = drv.prefix;
  parent_prefix_ = drv.parent_prefix;

  if (drv.ca) {
    process_continuous_assign_once(*drv.ca);
  } else {
    process_always_comb_once(*drv.comb);
  }

  prefix_ = saved_prefix;
  parent_prefix_ = saved_parent_prefix;
  resolving_wires_.erase(sym);
  return true;
}

Term SystemVerilogEncoder::lookup_symbol(const slang::ast::Symbol * sym)
{
  using namespace slang::ast;

  // Procedural for-loop counter: bound to a per-iteration BV
  // constant for the duration of the unrolling.
  auto lvt = loop_var_terms_.find(sym);
  if (lvt != loop_var_terms_.end()) {
    return lvt->second;
  }

  // If `sym` is a child instance's output-port internal, reconstruct
  // its value from the (one, in the common case; more than one for a
  // concatenation-target connection) segment(s) it was split across.
  // This may chase through multiple levels of instantiation (e.g. a
  // grandchild's output port connected straight through an
  // intermediate module's own output port, or one element of an
  // instance array wired to a slice of a parent-side bus), each
  // segment resolved all the way to its own non-aliased root.
  if (port_output_aliases_.count(sym)) {
    uint64_t width = sym->as<ValueSymbol>().getType().getBitWidth();
    auto pieces = resolve_output_alias_pieces(sym, 0, width - 1);
    std::sort(pieces.begin(),
              pieces.end(),
              [](const ResolvedAliasPiece & a, const ResolvedAliasPiece & b) {
                return a.rhs_lo > b.rhs_lo;
              });
    Term result;
    for (auto & piece : pieces) {
      Term t = lookup_symbol(piece.sym);
      Term piece_term = slice_bits(t, piece.target_lo, piece.target_hi);
      result =
          result ? solver_->make_term(Concat, result, piece_term) : piece_term;
    }
    return result;
  }

  // Wire being defined in the enclosing always_comb block: return
  // the partial accumulated term so that read-modify-write patterns
  // (e.g. `popcount = popcount + din[i];` inside an unrolled for
  // loop) see the previously-written value.
  auto pending_it = pending_comb_updates_.find(sym);
  if (pending_it != pending_comb_updates_.end()) {
    return pending_it->second;
  }

  auto it = symbol_to_term_.find(sym);
  if (it != symbol_to_term_.end()) {
    return it->second;
  }

  // Not resolved yet -- if `sym` is a wire whose driving continuous
  // assign / always_comb block simply hasn't been walked yet (e.g. it
  // appears later in program order than this read), process it now,
  // out of order, and retry.
  if (resolve_wire_on_demand(sym)) {
    auto pit = pending_comb_updates_.find(sym);
    if (pit != pending_comb_updates_.end()) return pit->second;
    auto sit = symbol_to_term_.find(sym);
    if (sit != symbol_to_term_.end()) return sit->second;
  }

  // Parameter / localparam: slang has already evaluated the value at
  // elaboration time.  Materialize a fresh BV constant from it so
  // references to `REQUESTERS`, `MAX_COUNT`, etc. fold to literals.
  if (sym->kind == SymbolKind::Parameter) {
    auto & param = sym->as<ParameterSymbol>();
    const auto & cv = param.getValue();
    if (!cv.isInteger()) {
      throw PonoException("SystemVerilogEncoder: non-integer parameter '"
                          + string(sym->name) + "'");
    }
    auto val = cv.integer();
    uint64_t width = param.getType().getBitWidth();
    if (width == 0) width = val.getBitWidth();
    if (width == 0) width = 32;
    Sort sort = solver_->make_sort(BV, width);
    val.setSigned(false);
    string val_str =
        val.toString(slang::LiteralBase::Decimal, /*includeBase=*/false);
    return solver_->make_term(val_str, sort, 10);
  }

  // Enum literal (`IDLE`, `REQ`, ...): slang has already evaluated its
  // value at elaboration time, exactly like a parameter.  Declaring an
  // enum-typed *variable* already worked (its type is integral, so
  // type_to_sort() succeeds); only referencing one of the enum's own
  // named values -- ordinary once the enum is declared -- was missing.
  if (sym->kind == SymbolKind::EnumValue) {
    auto & enum_val = sym->as<EnumValueSymbol>();
    const auto & cv = enum_val.getValue();
    if (!cv.isInteger()) {
      throw PonoException("SystemVerilogEncoder: non-integer enum value '"
                          + string(sym->name) + "'");
    }
    auto val = cv.integer();
    uint64_t width = enum_val.getType().getBitWidth();
    if (width == 0) width = val.getBitWidth();
    if (width == 0) width = 32;
    Sort sort = solver_->make_sort(BV, width);
    val.setSigned(false);
    string val_str =
        val.toString(slang::LiteralBase::Decimal, /*includeBase=*/false);
    return solver_->make_term(val_str, sort, 10);
  }

  throw PonoException("SystemVerilogEncoder: unknown symbol '"
                      + string(sym->name) + "'");
}

string SystemVerilogEncoder::make_name(const string & name) const
{
  if (prefix_.empty()) return name;
  return prefix_ + "." + name;
}

Term SystemVerilogEncoder::wire_seed_term(const slang::ast::Symbol * sym)
{
  auto it = symbol_to_term_.find(sym);
  if (it != symbol_to_term_.end()) return it->second;
  uint64_t width = sym->as<slang::ast::ValueSymbol>().getType().getBitWidth();
  if (width == 0) width = 1;
  Term iv = fts_.make_inputvar(make_name(string(sym->name)),
                               solver_->make_sort(BV, width));
  symbol_to_term_[sym] = iv;
  return iv;
}

Term SystemVerilogEncoder::resize_to(const Term & t,
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
    return solver_->make_term(ext_op, t);
  }
  // Truncate (extract lower bits).
  return solver_->make_term(Op(Extract, target_width - 1, 0), t);
}

Term SystemVerilogEncoder::replace_bits(const Term & base,
                                        const Term & slice,
                                        uint64_t lo,
                                        uint64_t hi)
{
  uint64_t base_w = base->get_sort()->get_width();
  if (lo == 0 && hi == base_w - 1) return slice;
  std::vector<Term> parts;
  if (hi + 1 < base_w) {
    parts.push_back(solver_->make_term(Op(Extract, base_w - 1, hi + 1), base));
  }
  parts.push_back(slice);
  if (lo > 0) {
    parts.push_back(solver_->make_term(Op(Extract, lo - 1, 0), base));
  }
  Term result = parts[0];
  for (size_t i = 1; i < parts.size(); ++i) {
    result = solver_->make_term(Concat, result, parts[i]);
  }
  return result;
}

Term SystemVerilogEncoder::replace_bits_dynamic(const Term & base,
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
  Sort base_sort = solver_->make_sort(BV, base_w);
  Term idx_ext = resize_to(idx, base_w, false);
  Term shift_amount = idx_ext;
  if (elem_w != 1) {
    Term elem_w_term = solver_->make_term(elem_w, base_sort);
    shift_amount = solver_->make_term(BVMul, idx_ext, elem_w_term);
  }
  Term elem_ones = solver_->make_term(
      BVNot, solver_->make_term(0, solver_->make_sort(BV, elem_w)));
  Term mask = solver_->make_term(
      BVShl, resize_to(elem_ones, base_w, false), shift_amount);
  Term slice_shifted =
      solver_->make_term(BVShl, resize_to(slice, base_w, false), shift_amount);
  Term cleared =
      solver_->make_term(BVAnd, base, solver_->make_term(BVNot, mask));
  return solver_->make_term(
      BVOr, cleared, solver_->make_term(BVAnd, slice_shifted, mask));
}

}  // namespace pono
