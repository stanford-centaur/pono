/*!
 * \file terms.cpp
 * \brief Low-level symbol/term-lookup helpers shared across the SV encoder.
 * \author Áron Ricardo Perez-Lopez
 * \date 2026
 * \copyright See the LICENSE file in the top-level source directory.
 *
 * Covers lookup_symbol()'s resolution chain (loop-variable bindings,
 * output-port-alias reconstruction via resolve_output_alias_pieces(),
 * in-progress always_comb partial values, on-demand wire resolution with
 * combinational-loop detection, and parameter/enum-literal materialization
 * from slang's elaborated constants), plus naming/lookup utilities
 * (make_name(), wire_seed_term()). The pure bit/type-manipulation helpers
 * (type_to_sort(), slice_bits(), resize_to(), replace_bits(),
 * replace_bits_dynamic()) live in bit_utils.h/.cpp instead, since they need
 * no encoder state at all.
 */
#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

#include "frontends/systemverilog/bit_utils.h"
#include "frontends/systemverilog/encoder.h"
#include "slang/ast/Symbol.h"
#include "slang/ast/symbols/ParameterSymbols.h"
#include "slang/ast/symbols/ValueSymbol.h"
#include "slang/ast/types/AllTypes.h"
#include "smt-switch/smt.h"
#include "utils/exceptions.h"

using namespace smt;
using namespace std;

namespace pono {

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
      Term piece_term =
          slice_bits(solver_, t, piece.target_lo, piece.target_hi);
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

}  // namespace pono
