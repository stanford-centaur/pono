/*!
 * \file symbol_table.cpp
 * \brief Symbol classification, term binding, and on-demand wire lookup.
 * \author Áron Ricardo Perez-Lopez
 * \date 2026
 * \copyright See the LICENSE file in the top-level source directory.
 *
 * See symbol_table.h for what this class covers and why the classification
 * maps are exposed as direct accessors.
 */
#include "frontends/systemverilog/symbol_table.h"

#include <algorithm>
#include <unordered_set>

#include "frontends/systemverilog/ast_helpers.h"
#include "frontends/systemverilog/bit_utils.h"
#include "slang/ast/Expression.h"
#include "slang/ast/SemanticFacts.h"
#include "slang/ast/Statement.h"
#include "slang/ast/Symbol.h"
#include "slang/ast/expressions/AssignmentExpressions.h"
#include "slang/ast/expressions/MiscExpressions.h"
#include "slang/ast/expressions/OperatorExpressions.h"
#include "slang/ast/statements/ConditionalStatements.h"
#include "slang/ast/statements/LoopStatements.h"
#include "slang/ast/statements/MiscStatements.h"
#include "slang/ast/symbols/BlockSymbols.h"
#include "slang/ast/symbols/InstanceSymbols.h"
#include "slang/ast/symbols/ParameterSymbols.h"
#include "slang/ast/symbols/PortSymbols.h"
#include "slang/ast/symbols/ValueSymbol.h"
#include "slang/ast/types/AllTypes.h"
#include "smt-switch/smt.h"
#include "utils/exceptions.h"

using namespace smt;
using namespace std;

namespace pono {

SymbolTable::SymbolTable(FunctionalTransitionSystem & fts,
                         const smt::SmtSolver & solver)
    : fts_(fts), solver_(solver)
{
}

// ============================================================================
// Classification (formerly prescan.cpp)
// ============================================================================

namespace {

// Collect blocking-assignment targets, separating full-width LHSes
// (wire candidates) from partial LHSes (which must be state vars so
// the assignment handler's add_constraint slice constraints are valid).
void collect_blocking_targets(
    const slang::ast::Statement & stmt,
    std::unordered_set<const slang::ast::Symbol *> & full,
    std::unordered_set<const slang::ast::Symbol *> & partial)
{
  using namespace slang::ast;

  switch (stmt.kind) {
    case StatementKind::ExpressionStatement: {
      auto & es = stmt.as<ExpressionStatement>();
      auto & expr = es.expr;
      if (expr.kind == ExpressionKind::Assignment) {
        auto & assign = expr.as<AssignmentExpression>();
        auto & lhs = assign.left();
        if (lhs.kind == ExpressionKind::NamedValue) {
          full.insert(&canonicalize_modport_port(
              lhs.as<NamedValueExpression>().symbol));
        } else if (lhs.kind == ExpressionKind::HierarchicalValue) {
          full.insert(&canonicalize_modport_port(
              lhs.as<HierarchicalValueExpression>().symbol));
        } else if (auto * base = find_lhs_base(lhs)) {
          partial.insert(base);
        }
      }
      break;
    }
    case StatementKind::Block: {
      auto & block = stmt.as<BlockStatement>();
      auto & body = block.body;
      if (body.kind == StatementKind::List) {
        auto & list = body.as<StatementList>();
        for (auto * s : list.list) collect_blocking_targets(*s, full, partial);
      } else {
        collect_blocking_targets(body, full, partial);
      }
      break;
    }
    case StatementKind::Conditional: {
      auto & cond = stmt.as<ConditionalStatement>();
      collect_blocking_targets(cond.ifTrue, full, partial);
      if (cond.ifFalse) collect_blocking_targets(*cond.ifFalse, full, partial);
      break;
    }
    case StatementKind::Case: {
      auto & cs = stmt.as<CaseStatement>();
      for (auto & item : cs.items)
        collect_blocking_targets(*item.stmt, full, partial);
      if (cs.defaultCase)
        collect_blocking_targets(*cs.defaultCase, full, partial);
      break;
    }
    case StatementKind::Timed: {
      auto & ts = stmt.as<TimedStatement>();
      collect_blocking_targets(ts.stmt, full, partial);
      break;
    }
    case StatementKind::ForLoop: {
      // Recurse into the body so writes inside the
      // (compile-time-unrolled) loop are seen during pre-scan.
      auto & loop = stmt.as<ForLoopStatement>();
      collect_blocking_targets(loop.body, full, partial);
      break;
    }
    case StatementKind::WhileLoop:
      collect_blocking_targets(
          stmt.as<WhileLoopStatement>().body, full, partial);
      break;
    case StatementKind::DoWhileLoop:
      collect_blocking_targets(
          stmt.as<DoWhileLoopStatement>().body, full, partial);
      break;
    case StatementKind::RepeatLoop:
      collect_blocking_targets(
          stmt.as<RepeatLoopStatement>().body, full, partial);
      break;
    case StatementKind::ForeachLoop:
      collect_blocking_targets(
          stmt.as<ForeachLoopStatement>().body, full, partial);
      break;
    default: break;
  }
}

}  // namespace

void SymbolTable::pre_scan_always_ff(const slang::ast::Statement & body)
{
  collect_nonblocking_targets(body, state_var_symbols_);
}

void SymbolTable::pre_scan_state_vars(
    const slang::ast::InstanceBodySymbol & body, std::string & prefix)
{
  using namespace slang::ast;
  walk_members(body, prefix, [&](const Symbol & member) {
    if (member.kind == SymbolKind::ProceduralBlock) {
      auto & proc = member.as<ProceduralBlockSymbol>();
      if (proc.procedureKind == ProceduralBlockKind::AlwaysFF
          || proc.procedureKind == ProceduralBlockKind::Always) {
        pre_scan_always_ff(proc.getBody());
      } else if (proc.procedureKind == ProceduralBlockKind::AlwaysLatch) {
        pre_scan_always_latch(proc.getBody());
      } else if (proc.procedureKind == ProceduralBlockKind::Initial) {
        if (auto * forever_body = as_forever_event_body(proc.getBody())) {
          pre_scan_always_ff(*forever_body);
        }
      }
    } else if (member.kind == SymbolKind::Instance) {
      pre_scan_state_vars(member.as<InstanceSymbol>().body, prefix);
    }
  });
}

void SymbolTable::pre_scan_always_comb(
    const slang::ast::Statement & body,
    const slang::ast::ProceduralBlockSymbol & proc,
    const string & prefix,
    const string & parent_prefix)
{
  // Collect blocking-assign LHS targets.  Bases written full-width
  // become wires (macro-substituted); bases written through bit /
  // range selects become state vars instead so the assignment
  // handler can use slice-equality add_constraint constraints (which
  // requires the term to be a state var, not an input var).  Mixed
  // full+partial writes also go to state vars to keep the splice
  // semantics correct.
  std::unordered_set<const slang::ast::Symbol *> full, partial;
  collect_blocking_targets(body, full, partial);
  for (auto * sym : full) {
    if (state_var_symbols_.count(sym)) continue;
    if (partial.count(sym)) {
      state_var_symbols_.insert(sym);
    } else {
      wire_symbols_.insert(sym);
      wire_drivers_[sym] = { nullptr, &proc, prefix, parent_prefix };
    }
  }
  for (auto * sym : partial) {
    if (!wire_symbols_.count(sym)) {
      state_var_symbols_.insert(sym);
    }
  }
}

void SymbolTable::pre_scan_always_latch(const slang::ast::Statement & body)
{
  std::unordered_set<const slang::ast::Symbol *> full, partial;
  collect_blocking_targets(body, full, partial);
  for (auto * sym : full) state_var_symbols_.insert(sym);
  for (auto * sym : partial) state_var_symbols_.insert(sym);
}

void SymbolTable::pre_scan_instance(const slang::ast::InstanceSymbol & inst,
                                    string & prefix)
{
  using namespace slang::ast;

  // Each output (or inout) port of the child is driven by the child's
  // logic.  In the parent's view the connected expression is usually a
  // simple NamedValue (e.g. `child c (.sum(y))`), but for an instance-
  // array element it's a constant-index slice of a parent-side bus
  // (e.g. `fifo_data_out[i]`, already resolved by slang), or a
  // concatenation of several parent-side signals (`.sum({hi, lo})`,
  // splitting the port's bits across each) -- either way, every
  // underlying base symbol becomes a wire, unless it's already known
  // to be a register (state var), which takes priority. Non-constant
  // indices are not yet supported.
  for (auto * pc : inst.getPortConnections()) {
    if (!pc) continue;
    if (pc->port.kind != SymbolKind::Port) continue;
    auto & port = pc->port.as<PortSymbol>();
    if (port.direction != ArgumentDirection::Out
        && port.direction != ArgumentDirection::InOut)
      continue;
    auto * conn_expr = pc->getExpression();
    if (!conn_expr) continue;
    // Slang wraps an output-port connection as an Assignment whose
    // left-hand side is the parent-side expression being driven.
    // Unwrap it; the child's logic effectively writes through to the
    // LHS.
    if (conn_expr->kind == ExpressionKind::Assignment) {
      conn_expr = &conn_expr->as<AssignmentExpression>().left();
    }
    auto mark_wire = [&](const Expression & target) {
      auto * parent_sym = find_lhs_base(target);
      if (!parent_sym) return;
      if (!state_var_symbols_.count(parent_sym)) {
        wire_symbols_.insert(parent_sym);
      }
    };
    if (conn_expr->kind == ExpressionKind::Concatenation) {
      for (auto * operand :
           conn_expr->as<ConcatenationExpression>().operands()) {
        mark_wire(*operand);
      }
    } else {
      mark_wire(*conn_expr);
    }
  }

  // Recurse into nested instances so any wires further down the
  // hierarchy are visible to declare_variables.
  walk_members(inst.body, prefix, [&](const Symbol & m) {
    if (m.kind == SymbolKind::Instance) {
      pre_scan_instance(m.as<InstanceSymbol>(), prefix);
    }
  });
}

// ============================================================================
// Lookup / binding (formerly terms.cpp's stateful half)
// ============================================================================

string SymbolTable::make_name(const string & prefix, const string & name) const
{
  if (prefix.empty()) return name;
  return prefix + "." + name;
}

vector<ResolvedAliasPiece> SymbolTable::resolve_output_alias_pieces(
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

bool SymbolTable::resolve_wire_on_demand(const slang::ast::Symbol * sym)
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

  if (drv.ca) {
    driver_resolver_->resolve_continuous_assign(
        *drv.ca, drv.prefix, drv.parent_prefix);
  } else {
    driver_resolver_->resolve_always_comb(
        *drv.comb, drv.prefix, drv.parent_prefix);
  }

  resolving_wires_.erase(sym);
  return true;
}

Term SymbolTable::lookup_symbol(const slang::ast::Symbol * sym)
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

Term SymbolTable::wire_seed_term(const slang::ast::Symbol * sym,
                                 const string & prefix)
{
  auto it = symbol_to_term_.find(sym);
  if (it != symbol_to_term_.end()) return it->second;
  uint64_t width = sym->as<slang::ast::ValueSymbol>().getType().getBitWidth();
  if (width == 0) width = 1;
  Term iv = fts_.make_inputvar(make_name(prefix, string(sym->name)),
                               solver_->make_sort(BV, width));
  symbol_to_term_[sym] = iv;
  return iv;
}

}  // namespace pono
