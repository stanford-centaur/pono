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

#include "frontends/systemverilog/bit_utils.h"
#include "slang/ast/ASTVisitor.h"
#include "slang/ast/Scope.h"
#include "slang/ast/TimingControl.h"
#include "slang/ast/expressions/AssignmentExpressions.h"
#include "slang/ast/expressions/MiscExpressions.h"
#include "slang/ast/expressions/OperatorExpressions.h"
#include "slang/ast/expressions/SelectExpressions.h"
#include "slang/ast/statements/ConditionalStatements.h"
#include "slang/ast/statements/LoopStatements.h"
#include "slang/ast/statements/MiscStatements.h"
#include "slang/ast/symbols/BlockSymbols.h"
#include "slang/ast/symbols/InstanceSymbols.h"
#include "slang/ast/symbols/MemberSymbols.h"
#include "slang/ast/symbols/VariableSymbols.h"
#include "slang/ast/types/AllTypes.h"
#include "slang/ast/types/Type.h"
#include "slang/numeric/SVInt.h"
#include "utils/exceptions.h"

using namespace std;

namespace pono {

namespace {

// Recurses into a concatenation-target LHS (`{carry, sum} <= ...`) so
// every operand's own base symbol is classified, mirroring the
// per-operand handling already used by begin_write() (the write-time
// counterpart) and pre_scan_instance()'s output-port aliasing --
// find_lhs_base() itself can't do this since a concatenation has more
// than one base symbol and its return type is a single Symbol*.
void insert_nonblocking_lhs_targets(
    const slang::ast::Expression & lhs,
    std::unordered_set<const slang::ast::Symbol *> & targets)
{
  using namespace slang::ast;
  if (lhs.kind == ExpressionKind::Concatenation) {
    for (auto * operand : lhs.as<ConcatenationExpression>().operands()) {
      insert_nonblocking_lhs_targets(*operand, targets);
    }
    return;
  }
  if (auto * base = find_lhs_base(lhs)) {
    targets.insert(base);
  }
}

}  // namespace

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
          // for partial writes (`arr[i] <= ...`) we classify the base,
          // and for a concatenation-target write (`{carry, sum} <=
          // ...`) every operand's own base, recursing to support
          // nested concatenations.
          insert_nonblocking_lhs_targets(assign.left(), targets);
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

std::optional<LValueDesc> resolve_lvalue(
    const slang::ast::Expression & lhs,
    slang::ast::EvalContext & ctx,
    const slang::ast::ElementSelectExpression ** array_elem)
{
  using namespace slang::ast;
  switch (lhs.kind) {
    case ExpressionKind::NamedValue: {
      auto * sym =
          &canonicalize_modport_port(lhs.as<NamedValueExpression>().symbol);
      uint64_t w = value_width(*lhs.type);
      if (w == 0) {
        throw PonoException("SystemVerilogEncoder: zero-width lvalue '"
                            + string(sym->name) + "'");
      }
      return LValueDesc{ sym, 0, w - 1, w };
    }
    case ExpressionKind::HierarchicalValue: {
      auto * sym = &canonicalize_modport_port(
          lhs.as<HierarchicalValueExpression>().symbol);
      uint64_t w = value_width(*lhs.type);
      if (w == 0) {
        throw PonoException("SystemVerilogEncoder: zero-width lvalue '"
                            + string(sym->name) + "'");
      }
      return LValueDesc{ sym, 0, w - 1, w };
    }
    case ExpressionKind::ElementSelect: {
      auto & sel = lhs.as<ElementSelectExpression>();
      // An unpacked-array element is not a bit range of its base. A
      // caller that can turn it into a Store says so by passing
      // array_elem and gets the element as a synthetic base; one that
      // cannot gets nullopt -- see the contract note in ast_helpers.h.
      if (sel.value().type->getCanonicalType().kind
          == SymbolKind::FixedSizeUnpackedArrayType) {
        if (!array_elem) return std::nullopt;
        *array_elem = &sel;
        uint64_t elem_w = value_width(*lhs.type);
        if (elem_w == 0) {
          throw PonoException(
              "SystemVerilogEncoder: zero-width unpacked-array element "
              "lvalue");
        }
        return LValueDesc{ nullptr, 0, elem_w - 1, elem_w };
      }
      auto inner = resolve_lvalue(sel.value(), ctx, array_elem);
      if (!inner) return std::nullopt;
      auto idx_cv = sel.selector().eval(ctx);
      if (!idx_cv.isInteger()) return std::nullopt;
      // Signed, since a declared range may run below zero.
      auto idx_opt = idx_cv.integer().as<int64_t>();
      if (!idx_opt) {
        throw PonoException(
            "SystemVerilogEncoder: invalid constant element-select index");
      }
      uint64_t elem_w = value_width(*lhs.type);
      if (elem_w == 0) {
        throw PonoException(
            "SystemVerilogEncoder: zero-width element-select lvalue");
      }
      // The declared index is the bit offset only for an [n:0]
      // range; anything else has to be converted first.
      uint64_t idx = 0;
      const Type & base_type = sel.value().type->getCanonicalType();
      if (base_type.kind == SymbolKind::PackedArrayType) {
        if (!packed_element_ordinal(
                base_type.as<PackedArrayType>(), *idx_opt, idx)) {
          throw PonoException(
              "SystemVerilogEncoder: element-select index out of bounds");
        }
      } else {
        if (*idx_opt < 0) {
          throw PonoException(
              "SystemVerilogEncoder: element-select index out of bounds");
        }
        idx = static_cast<uint64_t>(*idx_opt);
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
      auto inner = resolve_lvalue(sel.value(), ctx, array_elem);
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
      auto inner = resolve_lvalue(ma.value(), ctx, array_elem);
      if (!inner) return std::nullopt;
      auto & field = ma.member.as<FieldSymbol>();
      uint64_t w = value_width(field.getType());
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

void walk_members(const slang::ast::Scope & scope,
                  std::string & prefix,
                  const std::function<void(const slang::ast::Symbol &)> & fn)
{
  using namespace slang::ast;
  for (auto & m : scope.members()) {
    if (m.kind == SymbolKind::GenerateBlockArray) {
      // Generate-for: walk each instantiated entry, pushing a
      // bracket-indexed prefix so per-iteration variables get
      // unique hierarchical names like "<top>.ctr[0].count".
      auto & arr = m.as<GenerateBlockArraySymbol>();
      std::string saved_prefix = prefix;
      std::string arr_name = std::string(arr.name);
      if (arr_name.empty()) arr_name = arr.getExternalName();
      for (auto * entry : arr.entries) {
        if (!entry || entry->isUninstantiated) continue;
        std::string idx_str;
        if (entry->arrayIndex) {
          auto idx = *entry->arrayIndex;
          idx.setSigned(false);
          idx_str =
              idx.toString(slang::LiteralBase::Decimal, /*includeBase=*/false);
        } else {
          idx_str = std::to_string(entry->constructIndex);
        }
        prefix = saved_prefix + "." + arr_name + "[" + idx_str + "]";
        walk_members(*entry, prefix, fn);
      }
      prefix = saved_prefix;
    } else if (m.kind == SymbolKind::InstanceArray) {
      // Arrayed instantiation (`mod inst[N-1:0] (...)`): slang
      // creates one child Instance per element, each with its own
      // port connections already sliced to the correct bus range
      // (see AssignmentExpressions.cpp's use of InstanceSymbol::
      // arrayPath). Flatten them into the same dispatch as a plain
      // Instance, pushing a bracket-indexed prefix so each element's
      // state gets a unique hierarchical name -- the elements
      // themselves are unnamed (only the array is named).
      auto & arr = m.as<InstanceArraySymbol>();
      std::string saved_prefix = prefix;
      std::string arr_name = std::string(arr.name);
      for (size_t i = 0; i < arr.elements.size(); ++i) {
        auto * element = arr.elements[i];
        if (!element) continue;
        prefix = saved_prefix + "." + arr_name + "[" + std::to_string(i) + "]";
        fn(*element);
      }
      prefix = saved_prefix;
    } else if (m.kind == SymbolKind::GenerateBlock) {
      // Generate-if / generate-case: a single block scope.  Push
      // its name (or slang's synthesized "genblkN") as the suffix.
      auto & gb = m.as<GenerateBlockSymbol>();
      if (gb.isUninstantiated) continue;
      std::string saved_prefix = prefix;
      std::string block_name = std::string(gb.name);
      if (block_name.empty()) block_name = gb.getExternalName();
      prefix = saved_prefix + "." + block_name;
      walk_members(gb, prefix, fn);
      prefix = saved_prefix;
    } else {
      fn(m);
    }
  }
}

bool is_concurrent_assertion_only(const slang::ast::Statement & body)
{
  using namespace slang::ast;

  switch (body.kind) {
    case StatementKind::ConcurrentAssertion: return true;
    case StatementKind::Block:
      return is_concurrent_assertion_only(body.as<BlockStatement>().body);
    case StatementKind::List: {
      for (auto * s : body.as<StatementList>().list) {
        if (!is_concurrent_assertion_only(*s)) return false;
      }
      return true;
    }
    case StatementKind::Empty: return true;
    default: return false;
  }
}

bool is_edge_triggered(const slang::ast::Statement & body)
{
  using namespace slang::ast;

  // slang wraps a procedural block's body in a Block whose own body is
  // the event-controlled statement; unwrap one level to reach it.
  const Statement * s = &body;
  if (s->kind == StatementKind::Block) {
    const Statement * inner = &s->as<BlockStatement>().body;
    if (inner->kind == StatementKind::List) {
      auto & list = inner->as<StatementList>();
      if (list.list.size() != 1) return false;
      inner = list.list[0];
    }
    s = inner;
  }
  if (s->kind != StatementKind::Timed) return false;

  auto is_edge = [](const TimingControl & tc) {
    return tc.kind == TimingControlKind::SignalEvent
           && tc.as<SignalEventControl>().edge != EdgeKind::None;
  };
  const TimingControl & timing = s->as<TimedStatement>().timing;
  if (timing.kind == TimingControlKind::EventList) {
    for (auto * ev : timing.as<EventListControl>().events) {
      if (is_edge(*ev)) return true;
    }
    return false;
  }
  return is_edge(timing);
}

namespace {

// Records every block-local read by an expression that the walk has
// not yet seen assigned.
struct UnassignedReadVisitor
    : slang::ast::ASTVisitor<UnassignedReadVisitor,
                             slang::ast::VisitFlags::Expressions>
{
  const std::unordered_set<const slang::ast::Symbol *> & assigned;
  std::unordered_set<const slang::ast::Symbol *> & out;

  void handle(const slang::ast::NamedValueExpression & e)
  {
    const slang::ast::Symbol & sym = e.symbol;
    if (is_block_local(sym) && !assigned.count(&sym)) out.insert(&sym);
  }
};

void scan_reads(const slang::ast::Expression & expr,
                const std::unordered_set<const slang::ast::Symbol *> & assigned,
                std::unordered_set<const slang::ast::Symbol *> & out)
{
  UnassignedReadVisitor v{ {}, assigned, out };
  expr.visit(v);
}

// Threads a "definitely assigned so far" set through the block in
// program order. A branch contributes only what *every* arm assigns,
// and a loop body contributes nothing, since it may run zero times.
void scan_hold_locals(const slang::ast::Statement & stmt,
                      std::unordered_set<const slang::ast::Symbol *> & assigned,
                      std::unordered_set<const slang::ast::Symbol *> & out)
{
  using namespace slang::ast;

  auto branch = [&](const Statement & s) {
    auto copy = assigned;
    scan_hold_locals(s, copy, out);
    return copy;
  };

  switch (stmt.kind) {
    case StatementKind::ExpressionStatement: {
      auto & expr = stmt.as<ExpressionStatement>().expr;
      if (expr.kind == ExpressionKind::Assignment) {
        auto & assign = expr.as<AssignmentExpression>();
        scan_reads(assign.right(), assigned, out);
        const Expression & lhs = assign.left();
        if (lhs.kind == ExpressionKind::NamedValue) {
          const Symbol & sym = lhs.as<NamedValueExpression>().symbol;
          // Only a whole-variable write makes it definitely
          // assigned; a partial one reads the rest of the variable.
          // Tracked for every symbol, not just block-local ones, so
          // that collect_definitely_assigned() can see module-level
          // targets too -- which cannot disturb the hold-local
          // result, since only a block-local read is ever reported.
          assigned.insert(&sym);
        } else {
          scan_reads(lhs, assigned, out);
        }
      } else {
        scan_reads(expr, assigned, out);
      }
      break;
    }
    case StatementKind::VariableDeclaration: {
      auto & vds = stmt.as<VariableDeclStatement>();
      if (auto * init = vds.symbol.getInitializer()) {
        scan_reads(*init, assigned, out);
        assigned.insert(&vds.symbol);
      }
      break;
    }
    case StatementKind::Block: {
      auto & body = stmt.as<BlockStatement>().body;
      if (body.kind == StatementKind::List) {
        for (auto * s : body.as<StatementList>().list) {
          scan_hold_locals(*s, assigned, out);
        }
      } else {
        scan_hold_locals(body, assigned, out);
      }
      break;
    }
    case StatementKind::List:
      for (auto * s : stmt.as<StatementList>().list) {
        scan_hold_locals(*s, assigned, out);
      }
      break;
    case StatementKind::Conditional: {
      auto & cond = stmt.as<ConditionalStatement>();
      for (auto & c : cond.conditions) scan_reads(*c.expr, assigned, out);
      auto t_assigned = branch(cond.ifTrue);
      if (!cond.ifFalse) break;
      auto f_assigned = branch(*cond.ifFalse);
      for (auto * sym : t_assigned) {
        if (f_assigned.count(sym)) assigned.insert(sym);
      }
      break;
    }
    case StatementKind::Case: {
      auto & cs = stmt.as<CaseStatement>();
      scan_reads(cs.expr, assigned, out);
      std::unordered_set<const Symbol *> common;
      bool first = true;
      for (auto & item : cs.items) {
        for (auto * e : item.expressions) scan_reads(*e, assigned, out);
        auto arm = branch(*item.stmt);
        if (first) {
          common = arm;
          first = false;
        } else {
          for (auto it = common.begin(); it != common.end();) {
            it = arm.count(*it) ? std::next(it) : common.erase(it);
          }
        }
      }
      // Without a default arm some value matches nothing, so no arm's
      // assignments are guaranteed.
      if (!cs.defaultCase) break;
      auto def = branch(*cs.defaultCase);
      if (first) {
        common = def;
      } else {
        for (auto it = common.begin(); it != common.end();) {
          it = def.count(*it) ? std::next(it) : common.erase(it);
        }
      }
      for (auto * sym : common) assigned.insert(sym);
      break;
    }
    case StatementKind::Timed:
      scan_hold_locals(stmt.as<TimedStatement>().stmt, assigned, out);
      break;
    case StatementKind::ForLoop: {
      // The iteration variable is supplied by the unrolling, not by
      // any assignment in the body, so it counts as assigned there.
      auto & loop = stmt.as<ForLoopStatement>();
      auto inner = assigned;
      for (auto * lv : loop.loopVars) inner.insert(lv);
      scan_hold_locals(loop.body, inner, out);
      break;
    }
    case StatementKind::ForeachLoop: {
      auto & loop = stmt.as<ForeachLoopStatement>();
      auto inner = assigned;
      for (auto & dim : loop.loopDims) {
        if (dim.loopVar) inner.insert(dim.loopVar);
      }
      scan_hold_locals(loop.body, inner, out);
      break;
    }
    case StatementKind::WhileLoop:
      branch(stmt.as<WhileLoopStatement>().body);
      break;
    case StatementKind::DoWhileLoop:
      branch(stmt.as<DoWhileLoopStatement>().body);
      break;
    case StatementKind::RepeatLoop:
      branch(stmt.as<RepeatLoopStatement>().body);
      break;
    default: break;
  }
}

}  // namespace

void collect_hold_locals(const slang::ast::Statement & body,
                         std::unordered_set<const slang::ast::Symbol *> & out)
{
  std::unordered_set<const slang::ast::Symbol *> assigned;
  scan_hold_locals(body, assigned, out);
}

void collect_definitely_assigned(
    const slang::ast::Statement & body,
    std::unordered_set<const slang::ast::Symbol *> & out)
{
  // Same walk, kept for what it threads rather than what it reports:
  // `assigned` is the set that survives every branch.
  std::unordered_set<const slang::ast::Symbol *> reads;
  scan_hold_locals(body, out, reads);
}

bool is_block_local(const slang::ast::Symbol & sym)
{
  using namespace slang::ast;
  const Scope * scope = sym.getParentScope();
  if (!scope) return false;
  SymbolKind owner = scope->asSymbol().kind;
  // A subroutine's formals, return value and locals live for one call
  // exactly as a block's temporaries live for one execution, and are
  // bound the same way while its body is inlined.
  return owner == SymbolKind::StatementBlock || owner == SymbolKind::Subroutine;
}

}  // namespace pono
