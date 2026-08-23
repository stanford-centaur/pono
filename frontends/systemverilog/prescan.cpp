/*! \file prescan.cpp
 *  \brief SystemVerilogEncoder's pre-scan pass: classifies which base
 *         symbols become wires vs. state vars before any variable is
 *         declared.
 */
#include <unordered_set>

#include "frontends/systemverilog/ast_helpers.h"
#include "frontends/systemverilog/encoder.h"
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
#include "slang/ast/symbols/PortSymbols.h"
#include "utils/exceptions.h"

using namespace smt;
using namespace std;

namespace pono {

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

// Thin wrapper around collect_nonblocking_targets(); called from
// pre_scan_state_vars() for every always_ff/always block in the design.
void SystemVerilogEncoder::pre_scan_always_ff(
    const slang::ast::Statement & body)
{
  collect_nonblocking_targets(body, state_var_symbols_);
}

void SystemVerilogEncoder::pre_scan_state_vars(
    const slang::ast::InstanceBodySymbol & body)
{
  using namespace slang::ast;
  walk_members(body, [&](const Symbol & member) {
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
      pre_scan_state_vars(member.as<InstanceSymbol>().body);
    }
  });
}

void SystemVerilogEncoder::pre_scan_always_comb(
    const slang::ast::Statement & body,
    const slang::ast::ProceduralBlockSymbol & proc)
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
      wire_drivers_[sym] = { nullptr, &proc, prefix_, parent_prefix_ };
    }
  }
  for (auto * sym : partial) {
    if (!wire_symbols_.count(sym)) {
      state_var_symbols_.insert(sym);
    }
  }
}

void SystemVerilogEncoder::pre_scan_always_latch(
    const slang::ast::Statement & body)
{
  std::unordered_set<const slang::ast::Symbol *> full, partial;
  collect_blocking_targets(body, full, partial);
  for (auto * sym : full) state_var_symbols_.insert(sym);
  for (auto * sym : partial) state_var_symbols_.insert(sym);
}

void SystemVerilogEncoder::pre_scan_instance(
    const slang::ast::InstanceSymbol & inst)
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
  walk_members(inst.body, [&](const Symbol & m) {
    if (m.kind == SymbolKind::Instance) {
      pre_scan_instance(m.as<InstanceSymbol>());
    }
  });
}

}  // namespace pono
