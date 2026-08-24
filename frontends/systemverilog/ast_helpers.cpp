/*!
 * \file ast_helpers.cpp
 * \brief LHS-resolution, traversal, and control-flow helper implementations.
 * \author Áron Ricardo Perez-Lopez
 * \date 2026
 * \copyright See the LICENSE file in the top-level source directory.
 *
 * These helpers fall into three groups. Statement traversal
 * (for_each_stmt_in_block, collect_nonblocking_targets) walks
 * block/conditional/case/loop bodies to find non-blocking-assignment targets
 * during pre-scan and process_instance(). LHS/lvalue resolution
 * (canonicalize_modport_port, find_lhs_base, resolve_lvalue, LValueDesc) maps
 * an assignment's left-hand side through modport indirection down to a base
 * Symbol and constant bit range. Compile-time control flow (LoopControlSignal,
 * as_forever_event_body) lets process_statement() model break/continue/disable
 * across unrolled loops and recognize the legacy `initial forever @(...) body`
 * idiom as equivalent to an always block.
 */
#include "frontends/systemverilog/ast_helpers.h"

#include <string>
#include <utility>

#include "slang/ast/expressions/AssignmentExpressions.h"
#include "slang/ast/expressions/MiscExpressions.h"
#include "slang/ast/expressions/SelectExpressions.h"
#include "slang/ast/statements/ConditionalStatements.h"
#include "slang/ast/statements/LoopStatements.h"
#include "slang/ast/statements/MiscStatements.h"
#include "slang/ast/symbols/MemberSymbols.h"
#include "slang/ast/symbols/VariableSymbols.h"
#include "slang/ast/types/Type.h"
#include "slang/numeric/SVInt.h"
#include "utils/exceptions.h"

using namespace std;

namespace pono {

void collect_nonblocking_targets(
    const slang::ast::Statement & stmt,
    std::unordered_set<const slang::ast::Symbol *> & targets)
{
  using namespace slang::ast;

  switch (stmt.kind) {
    case StatementKind::ExpressionStatement: {
      auto & es = stmt.as<ExpressionStatement>();
      auto & expr = es.expr;
      if (expr.kind == ExpressionKind::Assignment) {
        auto & assign = expr.as<AssignmentExpression>();
        if (assign.isNonBlocking()) {
          // The LHS of a non-blocking assignment is a state variable;
          // for partial writes (`arr[i] <= ...`) we classify the base.
          if (auto * base = find_lhs_base(assign.left())) {
            targets.insert(base);
          }
        }
      }
      break;
    }
    case StatementKind::Block: {
      auto & block = stmt.as<BlockStatement>();
      for_each_stmt_in_block(block, [&](const Statement & s) {
        collect_nonblocking_targets(s, targets);
      });
      break;
    }
    case StatementKind::Conditional: {
      auto & cond = stmt.as<ConditionalStatement>();
      collect_nonblocking_targets(cond.ifTrue, targets);
      if (cond.ifFalse) {
        collect_nonblocking_targets(*cond.ifFalse, targets);
      }
      break;
    }
    case StatementKind::Case: {
      auto & cs = stmt.as<CaseStatement>();
      for (auto & item : cs.items) {
        collect_nonblocking_targets(*item.stmt, targets);
      }
      if (cs.defaultCase) {
        collect_nonblocking_targets(*cs.defaultCase, targets);
      }
      break;
    }
    case StatementKind::Timed: {
      auto & ts = stmt.as<TimedStatement>();
      collect_nonblocking_targets(ts.stmt, targets);
      break;
    }
    case StatementKind::ForLoop: {
      // Recurse into the body so NB-assigned registers inside the
      // (compile-time-unrolled) loop are seen during pre-scan.
      auto & loop = stmt.as<ForLoopStatement>();
      collect_nonblocking_targets(loop.body, targets);
      break;
    }
    case StatementKind::WhileLoop:
      collect_nonblocking_targets(stmt.as<WhileLoopStatement>().body, targets);
      break;
    case StatementKind::DoWhileLoop:
      collect_nonblocking_targets(stmt.as<DoWhileLoopStatement>().body,
                                  targets);
      break;
    case StatementKind::RepeatLoop:
      collect_nonblocking_targets(stmt.as<RepeatLoopStatement>().body, targets);
      break;
    case StatementKind::ForeachLoop:
      collect_nonblocking_targets(stmt.as<ForeachLoopStatement>().body,
                                  targets);
      break;
    default:
      // Other statement types: nothing to extract.
      break;
  }
}

const slang::ast::Symbol & canonicalize_modport_port(
    const slang::ast::Symbol & sym)
{
  using namespace slang::ast;
  if (sym.kind == SymbolKind::ModportPort) {
    auto & mp = sym.as<ModportPortSymbol>();
    if (mp.internalSymbol) return *mp.internalSymbol;
  }
  return sym;
}

const slang::ast::Symbol * find_lhs_base(const slang::ast::Expression & lhs)
{
  using namespace slang::ast;
  switch (lhs.kind) {
    case ExpressionKind::NamedValue:
      return &canonicalize_modport_port(lhs.as<NamedValueExpression>().symbol);
    case ExpressionKind::HierarchicalValue:
      return &canonicalize_modport_port(
          lhs.as<HierarchicalValueExpression>().symbol);
    case ExpressionKind::ElementSelect:
      return find_lhs_base(lhs.as<ElementSelectExpression>().value());
    case ExpressionKind::RangeSelect:
      return find_lhs_base(lhs.as<RangeSelectExpression>().value());
    case ExpressionKind::MemberAccess:
      return find_lhs_base(lhs.as<MemberAccessExpression>().value());
    default: return nullptr;
  }
}

std::optional<LValueDesc> resolve_lvalue(const slang::ast::Expression & lhs,
                                         slang::ast::EvalContext & ctx)
{
  using namespace slang::ast;
  switch (lhs.kind) {
    case ExpressionKind::NamedValue: {
      auto * sym =
          &canonicalize_modport_port(lhs.as<NamedValueExpression>().symbol);
      uint64_t w = lhs.type->getBitWidth();
      if (w == 0) {
        throw PonoException("SystemVerilogEncoder: zero-width lvalue '"
                            + string(sym->name) + "'");
      }
      return LValueDesc{ sym, 0, w - 1, w };
    }
    case ExpressionKind::HierarchicalValue: {
      auto * sym = &canonicalize_modport_port(
          lhs.as<HierarchicalValueExpression>().symbol);
      uint64_t w = lhs.type->getBitWidth();
      if (w == 0) {
        throw PonoException("SystemVerilogEncoder: zero-width lvalue '"
                            + string(sym->name) + "'");
      }
      return LValueDesc{ sym, 0, w - 1, w };
    }
    case ExpressionKind::ElementSelect: {
      auto & sel = lhs.as<ElementSelectExpression>();
      auto inner = resolve_lvalue(sel.value(), ctx);
      if (!inner) return std::nullopt;
      auto idx_cv = sel.selector().eval(ctx);
      if (!idx_cv.isInteger()) return std::nullopt;
      auto idx_opt = idx_cv.integer().as<uint64_t>();
      if (!idx_opt) {
        throw PonoException(
            "SystemVerilogEncoder: invalid constant element-select index");
      }
      uint64_t idx = *idx_opt;
      uint64_t elem_w = lhs.type->getBitWidth();
      if (elem_w == 0) {
        throw PonoException(
            "SystemVerilogEncoder: zero-width element-select lvalue");
      }
      uint64_t lo = inner->lo + idx * elem_w;
      uint64_t hi = lo + elem_w - 1;
      if (hi > inner->hi) {
        throw PonoException(
            "SystemVerilogEncoder: element-select index out of bounds");
      }
      return LValueDesc{ inner->base, lo, hi, inner->base_w };
    }
    case ExpressionKind::RangeSelect: {
      // Constant range-select write (`w[7:4] <= ...`): both bounds must
      // be compile-time constants, mirroring the read-side logic in
      // expr_to_term()'s RangeSelect case (which likewise requires
      // constant bounds and normalizes hi/lo regardless of whether the
      // source wrote `[hi:lo]`, `[base +: width]`, or `[base -: width]`
      // -- slang's `.left()`/`.right()` already reflect the resolved
      // bounds either way). Unlike ElementSelect, there is no dynamic-
      // range-select write fallback anywhere in this encoder, so a
      // non-constant bound throws immediately instead of silently
      // dropping the write.
      auto & sel = lhs.as<RangeSelectExpression>();
      auto inner = resolve_lvalue(sel.value(), ctx);
      if (!inner) return std::nullopt;
      auto & left_expr = sel.left();
      auto & right_expr = sel.right();
      if (!left_expr.getConstant() || !right_expr.getConstant()) {
        throw PonoException(
            "SystemVerilogEncoder: non-constant range select lvalue bounds");
      }
      auto hi_opt = left_expr.getConstant()->integer().as<uint64_t>();
      auto lo_opt = right_expr.getConstant()->integer().as<uint64_t>();
      if (!hi_opt || !lo_opt) {
        throw PonoException(
            "SystemVerilogEncoder: invalid range select lvalue bounds");
      }
      uint64_t local_hi = *hi_opt;
      uint64_t local_lo = *lo_opt;
      if (local_hi < local_lo) swap(local_hi, local_lo);
      uint64_t lo = inner->lo + local_lo;
      uint64_t hi = inner->lo + local_hi;
      if (hi > inner->hi) {
        throw PonoException(
            "SystemVerilogEncoder: range select lvalue out of bounds");
      }
      return LValueDesc{ inner->base, lo, hi, inner->base_w };
    }
    case ExpressionKind::MemberAccess: {
      // Packed-struct/union field write (`s.field <= ...`): narrow the
      // inner base's range by the field's own bitOffset, mirroring the
      // read-side logic in expr_to_term()'s MemberAccess case. Additive
      // offsets compose correctly for nested access (`s.a.x`), since
      // each FieldSymbol's bitOffset is relative to its own immediately
      // enclosing struct/union type.
      auto & ma = lhs.as<MemberAccessExpression>();
      if (ma.member.kind != SymbolKind::Field) {
        throw PonoException(
            "SystemVerilogEncoder: unsupported member access lvalue on "
            + std::string(ma.member.name));
      }
      auto inner = resolve_lvalue(ma.value(), ctx);
      if (!inner) return std::nullopt;
      auto & field = ma.member.as<FieldSymbol>();
      uint64_t w = field.getType().getBitWidth();
      if (w == 0) {
        throw PonoException("SystemVerilogEncoder: zero-width field lvalue '"
                            + string(field.name) + "'");
      }
      uint64_t lo = inner->lo + field.bitOffset;
      uint64_t hi = lo + w - 1;
      if (hi > inner->hi) {
        throw PonoException("SystemVerilogEncoder: field lvalue out of bounds");
      }
      return LValueDesc{ inner->base, lo, hi, inner->base_w };
    }
    default:
      throw PonoException(
          "SystemVerilogEncoder: unsupported lvalue "
          "expression kind "
          + to_string(static_cast<int>(lhs.kind)));
  }
}

// `initial forever @(...) body` is a legacy structural spelling of
// `always @(...) body` -- unlike a general `forever` (which has no
// static iteration bound at all and is a genuine architectural
// boundary of this encoder's compile-time-unrolling model), this
// specific shape runs its (timing-controlled) body exactly once per
// pono-cycle, exactly like an always_ff/always block does. Returns the
// forever loop's own (Timed) body if `stmt` matches this shape
// (allowing the single-statement Block wrapper slang gives an
// `initial` block's top-level statement), nullptr otherwise.
const slang::ast::Statement * as_forever_event_body(
    const slang::ast::Statement & stmt)
{
  using namespace slang::ast;
  const Statement * s = &stmt;
  if (s->kind == StatementKind::Block) {
    auto & block = s->as<BlockStatement>();
    auto & inner = block.body;
    if (inner.kind == StatementKind::List) {
      auto & list = inner.as<StatementList>();
      if (list.list.size() != 1) return nullptr;
      s = list.list[0];
    } else {
      s = &inner;
    }
  }
  if (s->kind != StatementKind::ForeverLoop) return nullptr;
  auto & forever_stmt = s->as<ForeverLoopStatement>();
  if (forever_stmt.body.kind != StatementKind::Timed) return nullptr;
  return &forever_stmt.body;
}

}  // namespace pono
