/*!
 * \file statement_encoder.cpp
 * \brief The process_statement() switch encoding SV procedural statements.
 * \author Áron Ricardo Perez-Lopez
 * \date 2026
 * \copyright See the LICENSE file in the top-level source directory.
 */
#include "frontends/systemverilog/statement_encoder.h"

#include <algorithm>
#include <functional>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "frontends/systemverilog/assertion_walker.h"
#include "frontends/systemverilog/ast_helpers.h"
#include "frontends/systemverilog/bit_utils.h"
#include "frontends/systemverilog/expr_encoder.h"
#include "frontends/systemverilog/symbol_table.h"
#include "slang/ast/EvalContext.h"
#include "slang/ast/Expression.h"
#include "slang/ast/Patterns.h"
#include "slang/ast/SemanticFacts.h"
#include "slang/ast/Statement.h"
#include "slang/ast/Symbol.h"
#include "slang/ast/expressions/AssignmentExpressions.h"
#include "slang/ast/expressions/CallExpression.h"
#include "slang/ast/expressions/MiscExpressions.h"
#include "slang/ast/expressions/OperatorExpressions.h"
#include "slang/ast/expressions/SelectExpressions.h"
#include "slang/ast/statements/ConditionalStatements.h"
#include "slang/ast/statements/LoopStatements.h"
#include "slang/ast/statements/MiscStatements.h"
#include "slang/ast/symbols/BlockSymbols.h"
#include "slang/ast/symbols/SubroutineSymbols.h"
#include "slang/ast/symbols/VariableSymbols.h"
#include "slang/ast/types/AllTypes.h"
#include "slang/ast/types/Type.h"
#include "slang/numeric/SVInt.h"
#include "smt-switch/smt.h"
#include "utils/exceptions.h"
#include "utils/logger.h"

using namespace smt;
using namespace std;

namespace pono {

StatementEncoder::StatementEncoder(SymbolTable & symbol_table,
                                   ExprEncoder & expr_encoder,
                                   AssertionWalker & assertion_walker,
                                   FunctionalTransitionSystem & fts,
                                   const smt::SmtSolver & solver)
    : symbol_table_(symbol_table),
      expr_encoder_(expr_encoder),
      assertion_walker_(assertion_walker),
      fts_(fts),
      solver_(solver)
{
}

void StatementEncoder::inline_subroutine_body_no_return(
    const slang::ast::Statement & body, const string & prefix)
{
  bool saved_in = in_subroutine_;
  const slang::ast::Symbol * saved_ret = current_return_var_;
  in_subroutine_ = true;
  current_return_var_ = nullptr;
  try {
    process_statement(body,
                      StmtContext::COMBINATIONAL,
                      solver_->make_term(true),
                      prefix,
                      nullptr);
  }
  catch (...) {
    in_subroutine_ = saved_in;
    current_return_var_ = saved_ret;
    throw;
  }
  in_subroutine_ = saved_in;
  current_return_var_ = saved_ret;
}

void StatementEncoder::inline_subroutine_body(
    const slang::ast::Statement & body,
    const slang::ast::Symbol & return_var,
    const string & prefix)
{
  const slang::ast::Symbol * saved = current_return_var_;
  bool saved_in = in_subroutine_;
  current_return_var_ = &return_var;
  in_subroutine_ = true;
  // Every write a body like this makes is to one of its own locals,
  // which the local-write path handles whatever the context, so the
  // context below never decides anything; a write to anything else
  // is refused there rather than encoded against the caller's.
  try {
    process_statement(body,
                      StmtContext::COMBINATIONAL,
                      solver_->make_term(true),
                      prefix,
                      nullptr);
  }
  catch (...) {
    current_return_var_ = saved;
    in_subroutine_ = saved_in;
    throw;
  }
  current_return_var_ = saved;
  in_subroutine_ = saved_in;
}

void StatementEncoder::process_dynamic_write(
    const slang::ast::Expression & base_expr,
    const slang::ast::Expression & index_expr,
    const slang::ast::Type & write_type,
    bool scale_by_width,
    int64_t pos_bias,
    const slang::ast::Expression & rhs_expr,
    StmtContext ctx,
    const Term & condition,
    const string & prefix)
{
  using namespace slang::ast;

  // The select sits either directly on a variable (`base[idx] = rhs`)
  // or on a statically-resolvable sub-range of one (`p[2][idx] = rhs`,
  // `s.field[idx] = rhs`), in which case the write lands at that
  // range's own offset within the base.
  const Symbol * sym = nullptr;
  uint64_t base_offset = 0;
  const Expression & inner = base_expr;
  if (inner.kind == ExpressionKind::NamedValue
      || inner.kind == ExpressionKind::HierarchicalValue) {
    sym = &canonicalize_modport_port(
        (inner.kind == ExpressionKind::NamedValue)
            ? inner.as<NamedValueExpression>().symbol
            : inner.as<HierarchicalValueExpression>().symbol);
  } else {
    auto inner_desc = resolve_lvalue(inner, expr_encoder_.eval_ctx());
    if (!inner_desc || !inner_desc->base) {
      // Only a *second* runtime index is left: the write position is
      // then a product of two unknowns that this splice cannot name.
      throw PonoException(
          "SystemVerilogEncoder: a dynamic-index write whose base is itself "
          "runtime-indexed is not supported");
    }
    sym = inner_desc->base;
    base_offset = inner_desc->lo;
  }
  bool aliased = symbol_table_.port_output_aliases().count(sym) > 0;
  uint64_t sym_w = value_width(sym->as<ValueSymbol>().getType());
  auto pieces = symbol_table_.resolve_output_alias_pieces(sym, 0, sym_w - 1);
  if (pieces.empty()) return;

  uint64_t elem_w = value_width(write_type);
  if (elem_w == 0) elem_w = 1;

  // An initial write needs no fixed slice after all. The splice
  // below composes onto the variable's own term, so the constraint
  // reads "r equals r with these bits replaced" -- a tautology
  // everywhere else and a definition exactly where the write lands,
  // which is what pins the selected bits while leaving the rest of
  // the initial value free.

  Term idx = expr_encoder_.expr_to_term(index_expr, prefix);
  Term rhs = expr_encoder_.expr_to_term(rhs_expr, prefix);
  rhs = resize_to(solver_, rhs, elem_w, write_type.isSigned());

  // Where the write lands within `sym`, wide enough that the
  // arithmetic below cannot wrap for any position `sym` has.
  uint64_t pos_w = idx->get_sort()->get_width();
  while ((uint64_t{ 1 } << pos_w) < sym_w + elem_w) ++pos_w;
  Sort pos_sort = solver_->make_sort(BV, pos_w);
  Term pos_in_sym = resize_to(solver_, idx, pos_w, /*is_signed=*/false);
  // A declared index is the bit position only for an `[n:0]` range.
  // The read side rebases the same way (see expr_encoder's dynamic
  // element select); without it a write into, say, `logic [7:4] r`
  // lands four bits too high and silently misses.
  {
    const slang::ast::Type & bt = base_expr.type->getCanonicalType();
    if (bt.kind == SymbolKind::PackedArrayType) {
      auto & range = bt.as<PackedArrayType>().range;
      if (range.left >= range.right) {
        if (range.lower() != 0) {
          pos_in_sym = solver_->make_term(
              BVSub,
              pos_in_sym,
              solver_->make_term(static_cast<int64_t>(range.lower()),
                                 pos_sort));
        }
      } else {
        pos_in_sym = solver_->make_term(
            BVSub,
            solver_->make_term(static_cast<int64_t>(range.upper()), pos_sort),
            pos_in_sym);
      }
    }
  }
  if (scale_by_width && elem_w != 1) {
    pos_in_sym = solver_->make_term(
        BVMul, pos_in_sym, solver_->make_term(elem_w, pos_sort));
  }
  if (pos_bias > 0) {
    pos_in_sym = solver_->make_term(
        BVAdd,
        pos_in_sym,
        solver_->make_term(static_cast<uint64_t>(pos_bias), pos_sort));
  } else if (pos_bias < 0) {
    pos_in_sym = solver_->make_term(
        BVSub,
        pos_in_sym,
        solver_->make_term(static_cast<uint64_t>(-pos_bias), pos_sort));
  }
  if (base_offset != 0) {
    pos_in_sym = solver_->make_term(
        BVAdd, pos_in_sym, solver_->make_term(base_offset, pos_sort));
  }

  // Splice into one alias target, at `position` bits into it, when
  // `guard` says the write lands there at all.
  auto commit_to = [&](const Symbol * target,
                       const Term & position,
                       const Term & guard) {
    bool wire_comb = ctx == StmtContext::COMBINATIONAL
                     && symbol_table_.wire_symbols().count(target);
    Term prev_base;
    Term state_term;
    if (ctx == StmtContext::NEXT_STATE) {
      auto sit = symbol_table_.symbol_to_term().find(target);
      if (sit == symbol_table_.symbol_to_term().end()) {
        throw PonoException("SystemVerilogEncoder: dynamic-index write to '"
                            + string(target->name)
                            + "', which has no declared term");
      }
      state_term = sit->second;
      auto pit = symbol_table_.pending_next_updates().find(state_term);
      prev_base = (pit != symbol_table_.pending_next_updates().end())
                      ? pit->second
                      : state_term;
    } else {
      // Combinational, wire or not: compose onto whatever the block
      // has written so far. The single constraint binding a non-wire
      // is emitted when the block ends.
      auto pit = symbol_table_.pending_comb_updates().find(target);
      if (pit != symbol_table_.pending_comb_updates().end()) {
        prev_base = pit->second;
      } else {
        auto sit = symbol_table_.symbol_to_term().find(target);
        if (sit != symbol_table_.symbol_to_term().end()) {
          prev_base = sit->second;
        }
      }
    }
    if (!prev_base) return;

    Term combined = replace_bits_at(solver_, prev_base, rhs, position, elem_w);
    if (guard) combined = solver_->make_term(Ite, guard, combined, prev_base);
    if (condition != solver_->make_term(true)) {
      combined = solver_->make_term(Ite, condition, combined, prev_base);
    }
    if (ctx == StmtContext::NEXT_STATE) {
      symbol_table_.pending_next_updates()[state_term] = combined;
    } else {
      if (wire_comb && aliased) {
        symbol_table_.pending_comb_aliased().insert(target);
      }
      symbol_table_.pending_comb_updates()[target] = combined;
    }
  };

  // The ordinary case: one alias piece covering the whole symbol, so
  // the position within the target is the position within `sym` and
  // every write lands there.
  uint64_t first_w = value_width(pieces[0].sym->as<ValueSymbol>().getType());
  if (pieces.size() == 1 && pieces[0].rhs_lo == 0 && pieces[0].target_lo == 0
      && pieces[0].target_hi + 1 == first_w) {
    commit_to(pieces[0].sym, pos_in_sym, Term());
    return;
  }

  // Otherwise the symbol is spread across several parent-side
  // signals (a concatenation-target connection, or one element of an
  // instance array wired to a slice of a shared bus). Which of them
  // the write reaches is only known at runtime, so every piece takes
  // a guarded splice and at most one of the guards can hold.
  for (auto & piece : pieces) {
    if (piece.rhs_hi - piece.rhs_lo != piece.target_hi - piece.target_lo) {
      throw PonoException(
          "SystemVerilogEncoder: an output-port alias segment for '"
          + string(sym->name) + "' does not preserve its width");
    }
    // Lands here iff the whole element sits inside this segment.
    Term lo_term = solver_->make_term(piece.rhs_lo, pos_sort);
    Term hi_term = solver_->make_term(piece.rhs_hi - (elem_w - 1), pos_sort);
    Term guard =
        solver_->make_term(And,
                           solver_->make_term(BVUge, pos_in_sym, lo_term),
                           solver_->make_term(BVUle, pos_in_sym, hi_term));
    // Rebase onto the target: same distance from the segment's start.
    Term position = solver_->make_term(BVSub, pos_in_sym, lo_term);
    if (piece.target_lo != 0) {
      position = solver_->make_term(
          BVAdd, position, solver_->make_term(piece.target_lo, pos_sort));
    }
    commit_to(piece.sym, position, guard);
  }
}

smt::Term StatementEncoder::constant_array_value(
    const slang::ConstantValue & cv,
    const slang::ast::FixedSizeUnpackedArrayType & arr,
    const smt::Sort & array_sort)
{
  using namespace slang::ast;
  if (!cv.isUnpacked()) return Term();

  UnpackedArrayInfo info = unpacked_array_info(solver_, arr);
  auto elements = cv.elements();
  if (elements.empty()) return Term();

  // An element of a multi-dimensional array is itself an array, so
  // the pattern nests as far as the dimensions do.
  auto element_term = [&](const slang::ConstantValue & v) -> Term {
    if (v.isUnpacked()) {
      return constant_array_value(
          v,
          arr.elementType.getCanonicalType().as<FixedSizeUnpackedArrayType>(),
          info.element_sort);
    }
    auto svint = v.integer();
    svint.setSigned(false);
    return solver_->make_term(
        svint.toString(slang::LiteralBase::Decimal, /*includeBase=*/false),
        info.element_sort,
        10);
  };

  // A pattern lists its values left to right, which is the order of
  // the *declared* indices -- so for a descending range the first
  // value belongs to the highest index, and position counts the
  // opposite way from the normalized index the array sort uses.
  bool descending = arr.range.left >= arr.range.right;
  auto normalized = [&](size_t i) {
    return descending ? elements.size() - 1 - i : i;
  };

  // Fill with the first element, then Store only the ones that
  // differ -- a uniform pattern, which is the common case, stays a
  // single constant array.
  Term first = element_term(elements[0]);
  if (!first) return Term();
  Term filled = solver_->make_term(first, array_sort);
  for (size_t i = 1; i < elements.size(); ++i) {
    if (elements[i] == elements[0]) continue;
    Term value = element_term(elements[i]);
    if (!value) return Term();
    Term idx = solver_->make_term(normalized(i),
                                  solver_->make_sort(BV, info.index_width));
    filled = solver_->make_term(Store, filled, idx, value);
  }
  return filled;
}

smt::Term StatementEncoder::array_from_pattern(
    const slang::ast::Expression & expr,
    const slang::ast::FixedSizeUnpackedArrayType & arr,
    const Term & seed,
    const string & prefix)
{
  using namespace slang::ast;
  if (expr.kind != ExpressionKind::SimpleAssignmentPattern
      && expr.kind != ExpressionKind::StructuredAssignmentPattern
      && expr.kind != ExpressionKind::ReplicatedAssignmentPattern) {
    return Term();
  }

  UnpackedArrayInfo info = unpacked_array_info(solver_, arr);
  // The three pattern kinds share a base, but it is abstract and so
  // not reachable through as<>(); each concrete one exposes the same
  // element list.
  std::span<const Expression * const> elements;
  switch (expr.kind) {
    case ExpressionKind::SimpleAssignmentPattern:
      elements = expr.as<SimpleAssignmentPatternExpression>().elements();
      break;
    case ExpressionKind::StructuredAssignmentPattern:
      elements = expr.as<StructuredAssignmentPatternExpression>().elements();
      break;
    default:
      elements = expr.as<ReplicatedAssignmentPatternExpression>().elements();
      break;
  }
  if (elements.empty() || elements.size() != info.depth) return Term();

  // Same ordering as the constant case: a pattern lists its values
  // in declared-index order, which for a descending range runs the
  // opposite way from the normalized index.
  bool descending = arr.range.left >= arr.range.right;
  auto normalized = [&](size_t i) {
    return descending ? elements.size() - 1 - i : i;
  };

  Sort idx_sort = solver_->make_sort(BV, info.index_width);
  Term filled = seed;
  for (size_t i = 0; i < elements.size(); ++i) {
    Term idx = solver_->make_term(normalized(i), idx_sort);
    Term value;
    if (info.element_sort->get_sort_kind() == ARRAY) {
      // An element of a multi-dimensional array is itself an array,
      // so a nested pattern recurses as far as the dimensions do --
      // seeded with the element it is replacing.
      value = array_from_pattern(
          *elements[i],
          arr.elementType.getCanonicalType().as<FixedSizeUnpackedArrayType>(),
          solver_->make_term(Select, seed, idx),
          prefix);
    } else {
      value = expr_encoder_.expr_to_term(*elements[i], prefix);
    }
    if (!value) return Term();
    filled = solver_->make_term(Store, filled, idx, value);
  }
  return filled;
}

smt::Term StatementEncoder::constant_array_term(
    const slang::ast::Expression & rhs_expr,
    const slang::ast::FixedSizeUnpackedArrayType & arr,
    const smt::Sort & array_sort)
{
  return constant_array_value(
      rhs_expr.eval(expr_encoder_.eval_ctx()), arr, array_sort);
}

smt::Term StatementEncoder::array_pending_value(const slang::ast::Symbol * sym,
                                                const Term & state_term,
                                                StmtContext ctx)
{
  // A clocked block accumulates by term, the others by symbol; both
  // start from the array's own term when nothing has written it yet.
  if (ctx == StmtContext::NEXT_STATE) {
    auto it = symbol_table_.pending_next_updates().find(state_term);
    return it == symbol_table_.pending_next_updates().end() ? state_term
                                                            : it->second;
  }
  auto it = symbol_table_.pending_comb_updates().find(sym);
  return it == symbol_table_.pending_comb_updates().end() ? state_term
                                                          : it->second;
}

void StatementEncoder::record_array_pending(const slang::ast::Symbol * sym,
                                            const Term & state_term,
                                            StmtContext ctx,
                                            const Term & value)
{
  if (ctx == StmtContext::NEXT_STATE) {
    symbol_table_.pending_next_updates()[state_term] = value;
  } else {
    symbol_table_.pending_comb_updates()[sym] = value;
  }
}

bool StatementEncoder::process_array_element_assign(
    const slang::ast::Expression & lhs_expr,
    const slang::ast::Expression & rhs_expr,
    StmtContext ctx,
    const Term & condition,
    const string & prefix,
    const Term & rhs_override)
{
  using namespace slang::ast;

  // Check the base symbol first: resolve_lvalue() throws on shapes
  // begin_write() handles itself (a concatenation target, say), so it
  // may only be consulted once the target is known to be an array.
  const Symbol * base = find_lhs_base(lhs_expr);
  if (!base) return false;
  auto * base_value = base->as_if<ValueSymbol>();
  if (!base_value
      || base_value->getType().getCanonicalType().kind
             != SymbolKind::FixedSizeUnpackedArrayType) {
    return false;
  }

  // resolve_lvalue() reports the element select it bottomed out at and
  // describes the written range *within* that element, so a plain
  // `mem[i]` and a narrower `mem[i][3:0]` / `mem[i].f` / `mem[i][2]`
  // all arrive here the same way.
  const ElementSelectExpression * elem_sel = nullptr;
  std::optional<LValueDesc> desc;
  // An element that is itself an array (`m[i] <= row`, where m has a
  // further dimension) is no bit range of anything, so resolve_lvalue()
  // cannot describe it; the whole element is simply replaced.
  bool elem_is_array = lhs_expr.kind == ExpressionKind::ElementSelect
                       && lhs_expr.type->getCanonicalType().kind
                              == SymbolKind::FixedSizeUnpackedArrayType
                       && lhs_expr.as<ElementSelectExpression>()
                                  .value()
                                  .type->getCanonicalType()
                                  .kind
                              == SymbolKind::FixedSizeUnpackedArrayType;
  if (elem_is_array) {
    elem_sel = &lhs_expr.as<ElementSelectExpression>();
  } else {
    desc = resolve_lvalue(lhs_expr, expr_encoder_.eval_ctx(), &elem_sel);
    if (!elem_sel) return false;
  }

  // What `elem_sel` indexes may itself be an element of an outer
  // dimension (`m[i][j]` indexes `m[i]`), so peel those too. Each one
  // becomes a Select on the way down and a Store on the way back out.
  const Expression * base_expr = &elem_sel->value();
  std::vector<const ElementSelectExpression *> outer;
  while (base_expr->kind == ExpressionKind::ElementSelect
         && base_expr->as<ElementSelectExpression>()
                    .value()
                    .type->getCanonicalType()
                    .kind
                == SymbolKind::FixedSizeUnpackedArrayType) {
    auto & sel = base_expr->as<ElementSelectExpression>();
    outer.push_back(&sel);
    base_expr = &sel.value();
  }
  std::reverse(outer.begin(), outer.end());

  if (base_expr->kind != ExpressionKind::NamedValue
      && base_expr->kind != ExpressionKind::HierarchicalValue) {
    throw PonoException(
        "SystemVerilogEncoder: an unpacked-array element write must select on "
        "a declared array directly");
  }
  const Symbol * sym = &canonicalize_modport_port(
      (base_expr->kind == ExpressionKind::NamedValue)
          ? base_expr->as<NamedValueExpression>().symbol
          : base_expr->as<HierarchicalValueExpression>().symbol);

  // The only shape resolve_lvalue() declines without throwing is a
  // runtime-variable bit position above the element (`mem[i][j]`),
  // which is a dynamic splice within the element rather than a fixed
  // range of it. Whatever sits between the element and that position
  // still resolves, and gives the offset to splice at.
  const ElementSelectExpression * dyn_bit = nullptr;
  uint64_t dyn_offset = 0;
  if (!desc && !elem_is_array) {
    if (lhs_expr.kind == ExpressionKind::ElementSelect) {
      auto & outer = lhs_expr.as<ElementSelectExpression>();
      const ElementSelectExpression * within = nullptr;
      auto inner =
          resolve_lvalue(outer.value(), expr_encoder_.eval_ctx(), &within);
      if (inner && within == elem_sel) {
        dyn_bit = &outer;
        dyn_offset = inner->lo;
        desc = inner;
      }
    }
    if (!dyn_bit) {
      throw PonoException(
          "SystemVerilogEncoder: writing a runtime-variable bit position "
          "inside an unpacked-array element ('"
          + string(sym->name) + "') is only supported directly on the "
          + "element");
    }
  }
  auto sit = symbol_table_.symbol_to_term().find(sym);
  if (sit == symbol_table_.symbol_to_term().end()) {
    throw PonoException("SystemVerilogEncoder: write to unpacked array '"
                        + string(sym->name) + "' has no declared term");
  }
  Term state_term = sit->second;
  Term whole = array_pending_value(sym, state_term, ctx);

  // Walk down to the innermost array, remembering each level so the
  // Stores can be rebuilt outwards afterwards.
  struct Level
  {
    Term array;     ///< the array this level indexes
    Term index;     ///< normalized index into it
    Term in_range;  ///< null when no index can miss
  };
  std::vector<Level> levels;
  Term prev_base = whole;
  for (auto * sel : outer) {
    UnpackedArrayInfo dim = unpacked_array_info(
        solver_,
        sel->value().type->getCanonicalType().as<FixedSizeUnpackedArrayType>());
    Level level;
    level.array = prev_base;
    level.index = normalize_array_index(
        solver_,
        expr_encoder_.expr_to_term(sel->selector(), prefix),
        dim,
        &level.in_range);
    levels.push_back(level);
    prev_base = solver_->make_term(Select, level.array, level.index);
  }

  // The innermost dimension is the one `elem_sel` indexes, which is
  // the base's own type only when there are no outer dimensions.
  UnpackedArrayInfo info =
      unpacked_array_info(solver_,
                          elem_sel->value()
                              .type->getCanonicalType()
                              .as<FixedSizeUnpackedArrayType>());
  Term in_range;
  Term idx = normalize_array_index(
      solver_,
      expr_encoder_.expr_to_term(elem_sel->selector(), prefix),
      info,
      &in_range);

  if (dyn_bit) {
    // `mem[i][j] = v`: splice at a position only known at runtime,
    // offset by wherever the enclosing field or range starts.
    uint64_t write_w = dyn_bit->type->getBitWidth();
    if (write_w == 0) write_w = 1;
    Term value = rhs_override ? rhs_override
                              : expr_encoder_.expr_to_term(rhs_expr, prefix);
    value = resize_to(solver_, value, write_w, rhs_expr.type->isSigned());
    Term new_elem = replace_bits_dynamic(
        solver_,
        solver_->make_term(Select, prev_base, idx),
        value,
        expr_encoder_.expr_to_term(dyn_bit->selector(), prefix),
        write_w,
        dyn_offset);
    Term spliced = solver_->make_term(Store, prev_base, idx, new_elem);
    if (in_range) {
      spliced = solver_->make_term(Ite, in_range, spliced, prev_base);
    }
    // Rebuild the enclosing dimensions around the innermost Store.
    for (size_t k = levels.size(); k-- > 0;) {
      Term outer_store =
          solver_->make_term(Store, levels[k].array, levels[k].index, spliced);
      spliced = levels[k].in_range
                    ? solver_->make_term(
                          Ite, levels[k].in_range, outer_store, levels[k].array)
                    : outer_store;
    }
    record_array_pending(
        sym,
        state_term,
        ctx,
        (condition == solver_->make_term(true))
            ? spliced
            : solver_->make_term(Ite, condition, spliced, whole));
    return true;
  }

  Term new_elem;
  if (elem_is_array) {
    Sort want = type_to_sort(solver_, *lhs_expr.type);
    // An assignment pattern (`m[i] <= '{default: 0}`) is a value of
    // the element's array sort, which expr_to_term() has no case for.
    if (!rhs_override) {
      new_elem = constant_array_term(
          rhs_expr,
          lhs_expr.type->getCanonicalType().as<FixedSizeUnpackedArrayType>(),
          want);
    }
    if (!new_elem) {
      new_elem = rhs_override ? rhs_override
                              : expr_encoder_.expr_to_term(rhs_expr, prefix);
    }
    if (new_elem->get_sort() != want) {
      throw PonoException(
          "SystemVerilogEncoder: assigning to the array element '"
          + string(sym->name) + "[...]' needs a value of sort "
          + want->to_string() + ", not " + new_elem->get_sort()->to_string());
    }
  }
  uint64_t range_w = elem_is_array ? 0 : desc->hi - desc->lo + 1;
  Term rhs;
  if (!elem_is_array) {
    rhs = rhs_override ? rhs_override
                       : expr_encoder_.expr_to_term(rhs_expr, prefix);
    rhs = resize_to(solver_, rhs, range_w, rhs_expr.type->isSigned());
  }
  if (!elem_is_array)
    new_elem = (range_w == desc->base_w)
                   ? rhs
                   : replace_bits(solver_,
                                  solver_->make_term(Select, prev_base, idx),
                                  rhs,
                                  desc->lo,
                                  desc->hi);

  Term combined = solver_->make_term(Store, prev_base, idx, new_elem);
  if (in_range) {
    // The LRM ignores a write outside the declared range; letting it
    // through would corrupt whichever real cell the truncated index
    // happens to land on.
    combined = solver_->make_term(Ite, in_range, combined, prev_base);
  }
  // Rebuild the enclosing dimensions around the innermost Store.
  for (size_t k = levels.size(); k-- > 0;) {
    Term outer_store =
        solver_->make_term(Store, levels[k].array, levels[k].index, combined);
    combined = levels[k].in_range
                   ? solver_->make_term(
                         Ite, levels[k].in_range, outer_store, levels[k].array)
                   : outer_store;
  }
  record_array_pending(
      sym,
      state_term,
      ctx,
      (condition == solver_->make_term(true))
          ? combined
          : solver_->make_term(Ite, condition, combined, whole));
  return true;
}

bool StatementEncoder::process_whole_array_assign(
    const slang::ast::Expression & lhs_expr,
    const slang::ast::Expression & rhs_expr,
    StmtContext ctx,
    const Term & condition,
    const string & prefix)
{
  using namespace slang::ast;

  if (lhs_expr.kind != ExpressionKind::NamedValue
      && lhs_expr.kind != ExpressionKind::HierarchicalValue) {
    return false;
  }
  const slang::ast::Type & lhs_type = lhs_expr.type->getCanonicalType();
  if (lhs_type.kind != SymbolKind::FixedSizeUnpackedArrayType) return false;

  const Symbol * sym = &canonicalize_modport_port(
      (lhs_expr.kind == ExpressionKind::NamedValue)
          ? lhs_expr.as<NamedValueExpression>().symbol
          : lhs_expr.as<HierarchicalValueExpression>().symbol);

  auto sit = symbol_table_.symbol_to_term().find(sym);
  if (sit == symbol_table_.symbol_to_term().end()) return false;
  Term state_term = sit->second;

  Term combined = constant_array_term(rhs_expr,
                                      lhs_type.as<FixedSizeUnpackedArrayType>(),
                                      state_term->get_sort());
  if (!combined) {
    // A pattern whose values are not all elaboration-time constants
    // still has one term per element, so it builds as stores.
    combined = array_from_pattern(rhs_expr,
                                  lhs_type.as<FixedSizeUnpackedArrayType>(),
                                  state_term,
                                  prefix);
  }
  if (!combined) {
    // Not constant and not a pattern, so it has to be another array
    // of the same shape (`b <= a`): copying one array term into
    // another needs no element-by-element expansion.
    combined = expr_encoder_.expr_to_term(rhs_expr, prefix);
    if (combined->get_sort() != state_term->get_sort()) {
      throw PonoException(
          "SystemVerilogEncoder: whole-array assignment to '"
          + string(sym->name)
          + "' needs either a constant pattern or an array of the same "
            "shape");
    }
  }

  Term prev_base = array_pending_value(sym, state_term, ctx);
  record_array_pending(
      sym,
      state_term,
      ctx,
      (condition == solver_->make_term(true))
          ? combined
          : solver_->make_term(Ite, condition, combined, prev_base));
  return true;
}

Term StatementEncoder::pattern_match(
    const slang::ast::Pattern & pat,
    const Term & value,
    const string & prefix,
    std::vector<std::pair<const slang::ast::Symbol *, Term>> & bindings)
{
  using namespace slang::ast;

  switch (pat.kind) {
    case PatternKind::Wildcard:
      // `.*` places no constraint at all.
      return solver_->make_term(true);

    case PatternKind::Variable: {
      // `.name` matches anything and names what it matched, so the
      // arm's statement can read it.
      auto & vp = pat.as<VariablePattern>();
      bindings.emplace_back(&vp.variable, value);
      return solver_->make_term(true);
    }

    case PatternKind::Constant: {
      auto & cp = pat.as<ConstantPattern>();
      Term lit = expr_encoder_.expr_to_term(cp.expr, prefix);
      lit = resize_to(solver_,
                      lit,
                      value->get_sort()->get_width(),
                      cp.expr.type->isSigned());
      return solver_->make_term(Equal, value, lit);
    }

    case PatternKind::Structure: {
      // Each field pattern applies to that field's own bits, found
      // the same way a packed-struct member read finds them.
      auto & sp = pat.as<StructurePattern>();
      Term all;
      for (auto & fp : sp.patterns) {
        uint64_t w = value_width(fp.field->getType());
        if (w == 0) {
          throw PonoException(
              "SystemVerilogEncoder: a pattern over field '"
              + std::string(fp.field->name)
              + "' cannot be matched, since the field has no bits");
        }
        uint64_t lo = fp.field->bitOffset;
        Term slice = slice_bits(solver_, value, lo, lo + w - 1);
        Term one = pattern_match(*fp.pattern, slice, prefix, bindings);
        all = all ? solver_->make_term(And, all, one) : one;
      }
      return all ? all : solver_->make_term(true);
    }

    case PatternKind::Tagged:
      // A tagged union's discriminant is the thing being matched
      // here, and this encoder has no union representation carrying
      // one -- see the unpacked-union exclusion.
      throw PonoException(
          "SystemVerilogEncoder: a `tagged` pattern is not supported, since "
          "a tagged union's discriminant is not modeled");

    default:
      throw PonoException("SystemVerilogEncoder: unsupported pattern kind "
                          + std::to_string(static_cast<int>(pat.kind))
                          + " in a `case ... matches`");
  }
}

void StatementEncoder::refresh_loop_var_term(
    const slang::ast::ValueSymbol & sym)
{
  auto * cur = expr_encoder_.eval_ctx().findLocal(&sym);
  if (!cur || !cur->isInteger()) {
    throw PonoException("SystemVerilogEncoder: local variable '"
                        + string(sym.name) + "' lost its constant value");
  }
  auto svint = cur->integer();
  uint64_t width = sym.getType().getBitWidth();
  if (width == 0) width = svint.getBitWidth();
  if (width == 0) width = 32;
  Sort sort = solver_->make_sort(BV, width);
  svint.setSigned(false);
  string val_str = svint.toString(slang::LiteralBase::Decimal, false);
  symbol_table_.loop_var_terms()[&sym] = solver_->make_term(val_str, sort, 10);
}

void StatementEncoder::process_statement(
    const slang::ast::Statement & stmt,
    StmtContext ctx,
    const Term & condition,
    const string & prefix,
    const slang::ast::Expression * default_disable_expr)
{
  using namespace slang::ast;

  switch (stmt.kind) {
    case StatementKind::ExpressionStatement: {
      auto & es = stmt.as<ExpressionStatement>();
      auto & expr = es.expr;

      // A write to a plain compile-time-unrolled local (a `for`/
      // `while`/`repeat`/`foreach` scratch variable, or any other
      // `VariableDeclaration` local) is neither a wire nor a state
      // variable, so the SMT-term machinery below can't handle it: in
      // NEXT_STATE context begin_write() finds no declared term for it
      // and silently returns no writes, while in COMBINATIONAL/INITIAL
      // context the write instead reaches commit_write()'s "has no
      // declared term" exception -- delegate the whole expression to
      // slang's own constant evaluator instead, exactly as
      // `ForLoopStatement`'s own step expressions already do, and
      // refresh the mirrored SMT constant so later condition
      // evaluation (a `while`/`do`-`while` test, a `for`-loop bound)
      // sees the new value.
      {
        const Expression * lval_expr = nullptr;
        if (expr.kind == ExpressionKind::Assignment) {
          lval_expr = &expr.as<AssignmentExpression>().left();
        } else if (expr.kind == ExpressionKind::UnaryOp) {
          auto & unop = expr.as<UnaryExpression>();
          if (unop.op == UnaryOperator::Preincrement
              || unop.op == UnaryOperator::Postincrement
              || unop.op == UnaryOperator::Predecrement
              || unop.op == UnaryOperator::Postdecrement) {
            lval_expr = &unop.operand();
          }
        }
        if (lval_expr) {
          if (auto * base = find_lhs_base(*lval_expr)) {
            auto & vsym = base->as<ValueSymbol>();
            bool bound_local =
                expr_encoder_.eval_ctx().findLocal(&vsym) != nullptr;
            if (bound_local && !expr.eval(expr_encoder_.eval_ctx()).bad()) {
              // Still a compile-time constant, so keep folding it:
              // loop bounds and unrolled conditions depend on that.
              refresh_loop_var_term(vsym);
              break;
            }
            // A local with a term of its own is one the pre-scan
            // found is read on a path that never writes it, so it
            // holds its value and is a register; it takes the
            // ordinary write path below.
            if (in_subroutine_ && !is_block_local(vsym)) {
              // Inlining happens wherever the call appeared, so a
              // write reaching out of the subroutine would land in
              // whatever context that was, under none of the
              // conditions guarding the call.
              throw PonoException(
                  "SystemVerilogEncoder: an inlined subroutine writes '"
                  + string(base->name) + "', which is not local to it");
            }
            if ((bound_local || is_block_local(vsym))
                && !symbol_table_.symbol_to_term().count(&vsym)) {
              // A procedural temporary holding a runtime value. It has
              // no term of its own to write into, so binding it to the
              // value it now carries is what makes later reads of it
              // mean that value.
              if (expr.kind != ExpressionKind::Assignment) {
                throw PonoException(
                    "SystemVerilogEncoder: '++'/'--' on the local variable '"
                    + string(base->name)
                    + "' needs a compile-time-constant value");
              }
              auto & assign = expr.as<AssignmentExpression>();
              if (assign.left().kind != ExpressionKind::NamedValue) {
                throw PonoException(
                    "SystemVerilogEncoder: a partial write to the local "
                    "variable '"
                    + string(base->name) + "' is not supported");
              }
              Term val = expr_encoder_.expr_to_term(assign.right(), prefix);
              val = resize_to(solver_,
                              val,
                              vsym.getType().getBitWidth(),
                              assign.right().type->isSigned());
              auto & bound = symbol_table_.loop_var_terms();
              auto prev = bound.find(&vsym);
              if (condition != solver_->make_term(true)) {
                if (prev == bound.end()) {
                  throw PonoException(
                      "SystemVerilogEncoder: the local variable '"
                      + string(base->name)
                      + "' is written only under a runtime condition, so it "
                        "has no value on the other path");
                }
                val = solver_->make_term(Ite, condition, val, prev->second);
              }
              // A constant binding would now be stale.
              expr_encoder_.eval_ctx().deleteLocal(&vsym);
              bound[&vsym] = val;
              break;
            }
          }
        }
      }

      // Resolves `lhs_expr` to its base symbol/bit-range/output-alias
      // and the "previous" full-base term a write should compose onto
      // -- shared by plain/compound assignment (whose RHS may read
      // this via an implicit LValueReference) and by `++`/`--` (which
      // reads it directly as "the current value" -- see slice_of()
      // below). Returns an empty vector if lhs_expr isn't a shape
      // resolve_lvalue() handles (e.g. a dynamic-index element
      // select), in which case the caller may fall back to
      // process_dynamic_element_assign() for ElementSelect LHSes.
      // Ordinarily returns exactly one piece; a concatenation-target
      // output-port alias (`.port({hi, lo})`) splits a single write
      // into one piece per operand, each tagged with which slice of
      // the write's overall rhs value (`rhs_lo`/`rhs_hi`, 0-indexed
      // from the write's own low bit) it covers.
      struct LValueWrite
      {
        const Symbol * sym;
        bool aliased;
        bool has_range;  // true if this piece doesn't cover all of `sym`
        uint64_t lo, hi, slice_w;
        Term prev_base;  // may be null -- see call sites
        bool wire_comb;
        Term state_term;  // only valid when ctx == NEXT_STATE
        uint64_t rhs_lo, rhs_hi;
      };
      std::function<std::vector<LValueWrite>(const Expression &)> begin_write =
          [&](const Expression & lhs_expr) -> std::vector<LValueWrite> {
        // Concatenation-target LHS (`{hi, lo} <= ...;`): unlike a
        // range-/element-select, this has more than one base symbol,
        // so it can't be represented as a single LValueDesc. Recurse
        // into each operand (MSB-first, matching the write's overall
        // rhs numbering) and rebase its own piece(s)' rhs_lo/rhs_hi
        // (relative to that operand's own width) into the whole
        // concatenation's numbering -- exactly mirroring how a
        // concatenation-target *port connection* is split into one
        // OutputAliasSegment per operand.
        // A streaming concatenation target (`{>>{hi, lo}} <= a;`) is
        // the same positional split over its stream expressions --
        // what differs is only which bits of the RHS reach them,
        // which the caller sorts out by un-re-ordering the RHS
        // first. A Streaming expression's own type is void, so its
        // width comes from the bitstream.
        if (lhs_expr.kind == ExpressionKind::Concatenation
            || lhs_expr.kind == ExpressionKind::Streaming) {
          bool streaming = lhs_expr.kind == ExpressionKind::Streaming;
          std::vector<const Expression *> operands;
          uint64_t total_w = 0;
          if (streaming) {
            auto & sc = lhs_expr.as<StreamingConcatenationExpression>();
            total_w = sc.getBitstreamWidth();
            for (auto & stream : sc.streams()) {
              if (stream.withExpr) {
                throw PonoException(
                    "SystemVerilogEncoder: a `with` range in a streaming "
                    "concatenation target is not supported");
              }
              operands.push_back(stream.operand);
            }
          } else {
            total_w = lhs_expr.type->getBitWidth();
            for (auto * operand :
                 lhs_expr.as<ConcatenationExpression>().operands()) {
              operands.push_back(operand);
            }
          }
          if (total_w == 0) return {};
          std::vector<LValueWrite> writes;
          uint64_t covered = 0;
          for (const Expression * operand : operands) {
            // Ahead of the width check, which would otherwise report
            // a non-integral target as a zero-width one.
            if (streaming && !operand->type->isIntegral()) {
              throw PonoException(
                  "SystemVerilogEncoder: only an integral target can be "
                  "unpacked from a stream; '"
                  + std::string(operand->type->toString()) + "' cannot");
            }
            uint64_t seg_w = value_width(*operand->type);
            if (seg_w == 0 || covered + seg_w > total_w) {
              throw PonoException(
                  "SystemVerilogEncoder: a concatenation-target operand of "
                  "width "
                  + std::to_string(seg_w) + " does not fit the target's "
                  + std::to_string(total_w) + " bits");
            }
            auto sub = begin_write(*operand);
            if (sub.empty()) {
              // Returning no writes here would discard the operands
              // that *did* resolve along with this one, silently
              // dropping the whole assignment.
              throw PonoException(
                  "SystemVerilogEncoder: a concatenation-target write with a "
                  "runtime-indexed operand is not supported");
            }
            uint64_t seg_lo = total_w - covered - seg_w;
            for (auto & w : sub) {
              w.rhs_lo += seg_lo;
              w.rhs_hi += seg_lo;
              writes.push_back(w);
            }
            covered += seg_w;
          }
          return writes;
        }

        auto desc = resolve_lvalue(lhs_expr, expr_encoder_.eval_ctx());
        if (!desc) return {};
        bool base_aliased =
            symbol_table_.port_output_aliases().count(desc->base) > 0;
        auto pieces = symbol_table_.resolve_output_alias_pieces(
            desc->base, desc->lo, desc->hi);

        std::vector<LValueWrite> writes;
        writes.reserve(pieces.size());
        for (auto & piece : pieces) {
          const Symbol * sym = piece.sym;
          uint64_t lo = piece.target_lo;
          uint64_t hi = piece.target_hi;
          uint64_t sym_w = value_width(sym->as<ValueSymbol>().getType());
          bool has_range = !(lo == 0 && hi + 1 == sym_w);

          LValueWrite w{ sym,    base_aliased, has_range,   lo,
                         hi,     hi - lo + 1,  Term(),      false,
                         Term(), piece.rhs_lo, piece.rhs_hi };
          w.wire_comb = ctx == StmtContext::COMBINATIONAL
                        && symbol_table_.wire_symbols().count(sym);
          if (w.wire_comb) {
            auto pit = symbol_table_.pending_comb_updates().find(sym);
            if (pit != symbol_table_.pending_comb_updates().end()) {
              w.prev_base = pit->second;
            }
          } else if (ctx == StmtContext::NEXT_STATE) {
            auto sit = symbol_table_.symbol_to_term().find(sym);
            if (sit == symbol_table_.symbol_to_term().end()) {
              // All-or-nothing: if any piece of a (possibly split)
              // write can't resolve a state term, don't commit a
              // partial write for the others either.
              return {};
            }
            w.state_term = sit->second;
            auto pit = symbol_table_.pending_next_updates().find(w.state_term);
            w.prev_base = (pit != symbol_table_.pending_next_updates().end())
                              ? pit->second
                              : w.state_term;
          } else {
            // COMBINATIONAL non-wire or INITIAL: prev_base is the
            // current (constant) value of the LHS used only for
            // self-reference (compound assignment, ++/--).
            auto sit = symbol_table_.symbol_to_term().find(sym);
            if (sit != symbol_table_.symbol_to_term().end()) {
              w.prev_base = sit->second;
            }
          }
          writes.push_back(w);
        }
        return writes;
      };

      // Slices [lo, hi] out of `base`, or returns a null Term if
      // `base` itself is null (no previous write to reference yet).
      auto slice_of = [&](const Term & base, uint64_t lo, uint64_t hi) -> Term {
        if (!base) return Term();
        uint64_t pw = base->get_sort()->get_width();
        if (lo == 0 && hi == pw - 1) return base;
        return solver_->make_term(Op(Extract, hi, lo), base);
      };

      // Commits `rhs` (already the correct slice_w-wide final value)
      // as the write described by `w` -- shared by plain/compound
      // assignment and by ++/--, which differ only in how `rhs` was
      // computed.
      auto commit_write = [&](const LValueWrite & w, const Term & rhs) {
        if (w.wire_comb) {
          if (w.aliased) symbol_table_.pending_comb_aliased().insert(w.sym);
          // Compose new full-base value from prev_base + slice rhs.
          // On the very first write to this wire within the block
          // there is no prev_base yet; treat it as covering the
          // whole symbol only if it actually does, checked against
          // the symbol's declared width (not the write's own slice
          // width, which is trivially equal to itself) -- otherwise
          // seed a fresh placeholder to splice into, e.g. a
          // `for (i) arr[i] = ...;` pattern writing one element per
          // iteration.
          Term combined;
          if (w.prev_base) {
            combined = replace_bits(solver_, w.prev_base, rhs, w.lo, w.hi);
          } else {
            uint64_t sym_w = value_width(w.sym->as<ValueSymbol>().getType());
            bool full_write = !w.has_range && w.lo == 0 && w.hi + 1 == sym_w;
            combined =
                full_write
                    ? rhs
                    : replace_bits(solver_,
                                   symbol_table_.wire_seed_term(w.sym, prefix),
                                   rhs,
                                   w.lo,
                                   w.hi);
          }
          if (condition == solver_->make_term(true)) {
            symbol_table_.pending_comb_updates()[w.sym] = combined;
          } else {
            Term def = w.prev_base ? w.prev_base : combined;
            symbol_table_.pending_comb_updates()[w.sym] =
                solver_->make_term(Ite, condition, combined, def);
          }
          return;
        }

        auto it = symbol_table_.symbol_to_term().find(w.sym);
        if (it == symbol_table_.symbol_to_term().end()) {
          // Every non-wire write target should already have a term by
          // the time any statement is processed: locals are
          // intercepted earlier in this function (see findLocal()
          // above), a wire goes through the w.wire_comb branch above
          // instead, and every remaining declared symbol gets a term
          // from declare_variables_internal()/process_port() before
          // process_assignments() ever runs. Throw rather than
          // silently drop the write if that invariant somehow doesn't
          // hold.
          throw PonoException(
              "SystemVerilogEncoder: write to '" + string(w.sym->name)
              + "' has no declared term (out-of-order or unsupported "
                "target)");
        }
        Term lhs_term = it->second;
        uint64_t base_w = lhs_term->get_sort()->get_width();
        bool full_write = (w.lo == 0 && w.hi == base_w - 1);

        switch (ctx) {
          case StmtContext::NEXT_STATE: {
            Term combined =
                full_write
                    ? rhs
                    : replace_bits(solver_, w.prev_base, rhs, w.lo, w.hi);
            Term update;
            if (condition == solver_->make_term(true)) {
              update = combined;
            } else {
              update =
                  solver_->make_term(Ite, condition, combined, w.prev_base);
            }
            symbol_table_.pending_next_updates()[w.state_term] = update;
            break;
          }
          case StmtContext::COMBINATIONAL: {
            // Non-wire LHS (e.g. an output-port reg, or a base that
            // is also written partially). Accumulate the whole-base
            // value the way the wire branch above does, rather than
            // constraining each write's own slice as it happens:
            // separate per-write constraints all bind the same term
            // at once, so `p = 0; p[3] = 1;` would assert both p == 0
            // and p[3] == 1 and make the design vacuous. One
            // constraint per symbol is emitted when the block ends
            // (process_always_comb()), and reads within the block see
            // the accumulated value through lookup_symbol().
            auto pit = symbol_table_.pending_comb_updates().find(w.sym);
            Term prev = pit != symbol_table_.pending_comb_updates().end()
                            ? pit->second
                            : lhs_term;
            Term combined =
                full_write ? rhs : replace_bits(solver_, prev, rhs, w.lo, w.hi);
            symbol_table_.pending_comb_updates()[w.sym] =
                (condition == solver_->make_term(true))
                    ? combined
                    : solver_->make_term(Ite, condition, combined, prev);
            break;
          }
          case StmtContext::INITIAL: {
            // Accumulated per symbol for the reason the combinational
            // case is: separate per-write constraints all bind the
            // same term at once, so `x = 0; x[3] = 1;` would assert
            // both x == 0 and x[3] == 1 and leave no satisfiable
            // initial state at all. process_initial() emits one
            // constraint per symbol once the block ends.
            auto pit = symbol_table_.pending_comb_updates().find(w.sym);
            Term prev = pit != symbol_table_.pending_comb_updates().end()
                            ? pit->second
                            : lhs_term;
            Term combined =
                full_write ? rhs : replace_bits(solver_, prev, rhs, w.lo, w.hi);
            symbol_table_.pending_comb_updates()[w.sym] =
                (condition == solver_->make_term(true))
                    ? combined
                    : solver_->make_term(Ite, condition, combined, prev);
            break;
          }
        }
      };

      // Reassembles "the current value of the whole (possibly split)
      // lvalue" by concatenating each piece's own current-value slice
      // MSB-first (highest rhs_hi first) -- needed both for an
      // LValueReference inside a compound-assignment RHS and for
      // `++`/`--`'s "read the current value" step. Returns a null Term
      // if any piece has no previous value to read yet.
      auto reassemble_current =
          [&](const std::vector<LValueWrite> & ws) -> Term {
        std::vector<const LValueWrite *> ordered;
        ordered.reserve(ws.size());
        for (auto & w : ws) ordered.push_back(&w);
        std::sort(ordered.begin(),
                  ordered.end(),
                  [](const LValueWrite * a, const LValueWrite * b) {
                    return a->rhs_lo > b->rhs_lo;
                  });
        Term result;
        for (auto * w : ordered) {
          Term piece = slice_of(w->prev_base, w->lo, w->hi);
          if (!piece) return Term();
          result = result ? solver_->make_term(Concat, result, piece) : piece;
        }
        return result;
      };

      if (expr.kind == ExpressionKind::Assignment) {
        auto & assign = expr.as<AssignmentExpression>();
        auto & lhs_expr = assign.left();
        auto & rhs_expr = assign.right();

        // In a clocked block a blocking write is visible to later
        // reads in the same block; record the target so lookup_symbol()
        // hands them the pending value. Done after the write below, so
        // this assignment's own RHS still reads the old value.
        auto note_blocking_write = [&] {
          if (ctx != StmtContext::NEXT_STATE || assign.isNonBlocking()) return;
          if (auto * base = find_lhs_base(lhs_expr)) {
            symbol_table_.blocking_next_written().insert(base);
          }
        };

        // A whole-array assignment (`mem <= '0`) targets neither a bit
        // range nor a single element, so LValueDesc cannot describe it
        // and it is handled before begin_write().
        if (process_whole_array_assign(
                lhs_expr, rhs_expr, ctx, condition, prefix)) {
          note_blocking_write();
          break;
        }
        // Likewise for an unpacked-array element, or a bit range
        // inside one: a Store, never a commit_write().
        if (process_array_element_assign(
                lhs_expr, rhs_expr, ctx, condition, prefix)) {
          note_blocking_write();
          break;
        }

        auto writes = begin_write(lhs_expr);
        if (writes.empty()) {
          // resolve_lvalue() only handles constant-index selects; a
          // runtime-variable index (`arr[idx] = rhs`) needs a
          // dynamic-position splice instead of a static bit range.
          if (lhs_expr.kind == ExpressionKind::ElementSelect) {
            auto & sel = lhs_expr.as<ElementSelectExpression>();
            process_dynamic_write(sel.value(),
                                  sel.selector(),
                                  *sel.type,
                                  /*scale_by_width=*/true,
                                  /*pos_bias=*/0,
                                  rhs_expr,
                                  ctx,
                                  condition,
                                  prefix);
            note_blocking_write();
          } else if (lhs_expr.kind == ExpressionKind::RangeSelect) {
            // `r[i +: w]` and `r[i -: w]`: a fixed-width window at a
            // runtime position, which is the same splice an element
            // select needs, only already counted in bits. `-:` names
            // the top of its window, so the position is that many
            // bits lower.
            auto & rs = lhs_expr.as<RangeSelectExpression>();
            uint64_t w = value_width(*rs.type);
            bool down =
                rs.getSelectionKind() == RangeSelectionKind::IndexedDown;
            if (rs.getSelectionKind() == RangeSelectionKind::Simple) {
              throw PonoException(
                  "SystemVerilogEncoder: a range-select write with "
                  "non-constant `[hi:lo]` bounds is not supported; `+:` or "
                  "`-:` names a fixed width and is");
            }
            process_dynamic_write(rs.value(),
                                  rs.left(),
                                  *rs.type,
                                  /*scale_by_width=*/false,
                                  down ? -static_cast<int64_t>(w - 1) : 0,
                                  rhs_expr,
                                  ctx,
                                  condition,
                                  prefix);
            note_blocking_write();
          }
          break;
        }

        // Stash the slice value for any LValueReference inside rhs.
        Term saved_lvalue =
            expr_encoder_.set_current_lvalue_term(reassemble_current(writes));

        Term rhs = expr_encoder_.expr_to_term(rhs_expr, prefix);
        expr_encoder_.set_current_lvalue_term(saved_lvalue);

        uint64_t total_w = 0;
        for (auto & w : writes) total_w = std::max(total_w, w.rhs_hi + 1);
        if (lhs_expr.kind == ExpressionKind::Streaming) {
          // A stream is consumed from its most significant end, the
          // opposite of the truncation every other target wants, and
          // a source with too few bits is an error rather than
          // something to pad (LRM 11.4.14.3 -- slang rejects it
          // before this point). What is left is the stream the `<<`
          // re-ordering produced, so run that backwards to recover
          // the bits the targets are cut from.
          uint64_t rhs_w = rhs->get_sort()->get_width();
          if (rhs_w < total_w) {
            throw PonoException(
                "SystemVerilogEncoder: a streaming-concatenation target "
                "needs " + std::to_string(total_w) + " bits but the source "
                "supplies only " + std::to_string(rhs_w));
          }
          rhs = slice_bits(solver_, rhs, rhs_w - total_w, rhs_w - 1);
          rhs = stream_unreorder(
              solver_,
              rhs,
              lhs_expr.as<StreamingConcatenationExpression>().getSliceSize());
        } else {
          // A concatenation-target LHS is always unsigned per the LRM
          // (positional bit-splicing, not a numeric value); otherwise
          // use the RHS expression's own signedness.
          bool rhs_signed = lhs_expr.kind != ExpressionKind::Concatenation
                            && rhs_expr.type->isSigned();
          rhs = resize_to(solver_, rhs, total_w, rhs_signed);
        }
        for (auto & w : writes) {
          commit_write(w, slice_of(rhs, w.rhs_lo, w.rhs_hi));
        }
        note_blocking_write();
      } else if (expr.kind == ExpressionKind::UnaryOp) {
        // `i++`/`--i`/etc. as a standalone statement (distinct from
        // the same operators used as a `for`-loop step expression,
        // which slang's own constant evaluator already handles
        // separately via ForLoopStatement's step evaluation). Per the
        // LRM these are equivalent to `i = i +/- 1`; reuse the exact
        // same lvalue-resolution/commit machinery as plain assignment
        // above, reading the current value directly rather than
        // evaluating an RHS expression.
        auto & unop = expr.as<UnaryExpression>();
        if (unop.op == UnaryOperator::Preincrement
            || unop.op == UnaryOperator::Postincrement
            || unop.op == UnaryOperator::Predecrement
            || unop.op == UnaryOperator::Postdecrement) {
          auto writes = begin_write(unop.operand());
          if (writes.empty()) {
            const Symbol * base = find_lhs_base(unop.operand());
            auto * base_value = base ? base->as_if<ValueSymbol>() : nullptr;
            if (base_value
                && base_value->getType().getCanonicalType().kind
                       == SymbolKind::FixedSizeUnpackedArrayType) {
              // An array element is not a bit range of its base, so
              // begin_write() never describes one. Read it, step it,
              // and hand the result to the Store path.
              bool inc = unop.op == UnaryOperator::Preincrement
                         || unop.op == UnaryOperator::Postincrement;
              Term cur = expr_encoder_.expr_to_term(unop.operand(), prefix);
              Term stepped =
                  solver_->make_term(inc ? BVAdd : BVSub,
                                     cur,
                                     solver_->make_term(1, cur->get_sort()));
              if (process_array_element_assign(unop.operand(),
                                               unop.operand(),
                                               ctx,
                                               condition,
                                               prefix,
                                               stepped)) {
                break;
              }
            }
            logger.log(1,
                       "SystemVerilogEncoder: skipping unsupported ++/-- "
                       "operand shape");
            break;
          }
          Term cur = reassemble_current(writes);
          if (!cur) {
            throw PonoException(
                "SystemVerilogEncoder: '++'/'--' has no previous value to "
                "read for '"
                + std::string(writes[0].sym->name) + "'");
          }
          bool is_inc = unop.op == UnaryOperator::Preincrement
                        || unop.op == UnaryOperator::Postincrement;
          Term one = solver_->make_term(1, cur->get_sort());
          Term new_val = solver_->make_term(is_inc ? BVAdd : BVSub, cur, one);
          for (auto & w : writes) {
            commit_write(w, slice_of(new_val, w.rhs_lo, w.rhs_hi));
          }
        }
      } else if (expr.kind == ExpressionKind::Call) {
        // A bare call statement (`task_or_func(args);`, no assignment).
        // A system call (`$display`, `$finish`, `$fatal`, file I/O,
        // assertion/coverage control, etc.) has no synthesis meaning and
        // no effect on any modeled state, so it's safe to skip like the
        // other simulation-only constructs above (final blocks,
        // force/release, specify blocks). A user-defined task or
        // function has a body this encoder never inlines as a
        // statement, so any side effect it has on design state (writes
        // to variables read elsewhere) would be silently lost -- throw
        // instead of risking an unsound model.
        auto & call = expr.as<CallExpression>();
        if (call.isSystemCall()) {
          logger.log(1,
                     "SystemVerilogEncoder: skipping simulation-only system "
                     "call '{}' used as a statement",
                     call.getSubroutineName());
        } else {
          // A task's whole purpose is what it writes back, so inline
          // the body and then perform those writes here, in the
          // caller's own context and under the condition guarding the
          // call.
          auto * const * sub_ptr =
              std::get_if<const SubroutineSymbol *>(&call.subroutine);
          if (!sub_ptr || !*sub_ptr) {
            throw PonoException("SystemVerilogEncoder: unsupported call to '"
                                + std::string(call.getSubroutineName())
                                + "' used as a statement");
          }
          const SubroutineSymbol & sub = **sub_ptr;
          const string sub_name(sub.name);
          auto reject = [&](const string & why) {
            throw PonoException(
                "SystemVerilogEncoder: cannot inline the call to '" + sub_name
                + "': " + why);
          };
          if (inlining_tasks_.count(&sub)) reject("it is recursive");
          auto formals = sub.getArguments();
          if (formals.size() != call.arguments().size()) {
            reject(
                "it is called with a different number of arguments than "
                "it declares");
          }
          for (auto * formal : formals) {
            if (formal->direction == ArgumentDirection::Ref) {
              reject("argument '" + string(formal->name)
                     + "' is passed by reference");
            }
          }

          // slang hands an output/inout actual over as the copy-out
          // assignment itself, whose left side is the caller's own
          // variable -- both the read and the write-back want that.
          auto actual_of = [&](size_t k) -> const Expression & {
            const Expression * a = call.arguments()[k];
            if (a->kind == ExpressionKind::Assignment) {
              return a->as<AssignmentExpression>().left();
            }
            return *a;
          };

          auto & bound = symbol_table_.loop_var_terms();
          std::vector<std::pair<const Symbol *, Term>> saved_bindings;
          auto remember = [&](const Symbol * sym) {
            auto it = bound.find(sym);
            saved_bindings.emplace_back(
                sym, it == bound.end() ? Term() : it->second);
          };
          auto restore_bindings = [&]() {
            for (auto & entry : saved_bindings) {
              if (entry.second) {
                bound[entry.first] = entry.second;
              } else {
                bound.erase(entry.first);
              }
            }
          };

          // Read every actual first, in the caller's scope.
          std::vector<Term> incoming(formals.size());
          for (size_t k = 0; k < formals.size(); ++k) {
            if (formals[k]->direction == ArgumentDirection::Out) continue;
            uint64_t w = formals[k]->getType().getBitWidth();
            if (w == 0) {
              reject("argument '" + string(formals[k]->name)
                     + "' has no width");
            }
            const Expression & actual = actual_of(k);
            incoming[k] = resize_to(solver_,
                                    expr_encoder_.expr_to_term(actual, prefix),
                                    w,
                                    actual.type->isSigned());
          }
          for (size_t k = 0; k < formals.size(); ++k) {
            remember(formals[k]);
            if (incoming[k]) {
              bound[formals[k]] = incoming[k];
            } else {
              bound.erase(formals[k]);
            }
          }

          inlining_tasks_.insert(&sub);
          try {
            inline_subroutine_body_no_return(sub.getBody(), prefix);
          }
          catch (...) {
            inlining_tasks_.erase(&sub);
            restore_bindings();
            throw;
          }
          inlining_tasks_.erase(&sub);

          // Deliver each written-back argument to the caller's own
          // variable, collecting the values before any write so one
          // argument cannot be read back through another.
          std::vector<std::pair<size_t, Term>> outgoing;
          for (size_t k = 0; k < formals.size(); ++k) {
            if (formals[k]->direction == ArgumentDirection::In) continue;
            auto it = bound.find(formals[k]);
            if (it == bound.end()) {
              restore_bindings();
              reject("no path through it assigns argument '"
                     + string(formals[k]->name) + "'");
            }
            outgoing.emplace_back(k, it->second);
          }
          restore_bindings();

          for (auto & out : outgoing) {
            auto writes = begin_write(actual_of(out.first));
            if (writes.empty()) {
              reject("argument '" + string(formals[out.first]->name)
                     + "' is written back to something this encoder cannot "
                       "assign to");
            }
            uint64_t total_w = 0;
            for (auto & w : writes) total_w = std::max(total_w, w.rhs_hi + 1);
            Term value = resize_to(solver_,
                                   out.second,
                                   total_w,
                                   formals[out.first]->getType().isSigned());
            for (auto & w : writes) {
              commit_write(w, slice_of(value, w.rhs_lo, w.rhs_hi));
            }
          }
          break;
        }
      } else {
        throw PonoException(
            "SystemVerilogEncoder: unsupported expression statement kind "
            + std::to_string(static_cast<int>(expr.kind)));
      }
      break;
    }

    case StatementKind::List: {
      // A bare sequence of statements, with no block around it to
      // name or to absorb a `disable`. A subroutine body arrives this
      // way, so without this the body walks to nothing at all.
      for (auto * s : stmt.as<StatementList>().list) {
        process_statement(*s, ctx, condition, prefix, default_disable_expr);
      }
      break;
    }

    case StatementKind::Block: {
      auto & block = stmt.as<BlockStatement>();
      try {
        for_each_stmt_in_block(block, [&](const Statement & s) {
          process_statement(s, ctx, condition, prefix, default_disable_expr);
        });
      }
      catch (const LoopControlSignal & sig) {
        // Absorb a `disable <this block's name>;` reached from inside
        // (stopping the rest of this block); anything else -- a
        // Break/Continue meant for an enclosing ForLoop, or a Disable
        // targeting a different (typically outer) named block --
        // keeps propagating.
        if (sig.kind == LoopControlSignal::Disable && block.blockSymbol
            && sig.disable_target == block.blockSymbol) {
          break;
        }
        throw;
      }
      break;
    }

    case StatementKind::Conditional: {
      auto & cond_stmt = stmt.as<ConditionalStatement>();

      // If the condition is a compile-time constant (e.g. it only
      // references already-unrolled `for`-loop counters), branch on
      // it directly in C++ instead of building a symbolic guard for
      // both arms -- this is what lets break/continue/disable, which
      // can only be modeled as C++-level control flow
      // (LoopControlSignal), propagate correctly out of whichever
      // branch is actually taken. Skipped for the (rare) pattern-match
      // `if` form (multiple conditions); falls back to the general
      // symbolic-guard path below for any condition that isn't
      // const-evaluable (e.g. depends on a runtime signal).
      if (cond_stmt.conditions.size() == 1) {
        auto const_cv =
            cond_stmt.conditions[0].expr->eval(expr_encoder_.eval_ctx());
        if (!const_cv.bad()) {
          if (const_cv.isTrue()) {
            process_statement(
                cond_stmt.ifTrue, ctx, condition, prefix, default_disable_expr);
          } else if (cond_stmt.ifFalse) {
            process_statement(*cond_stmt.ifFalse,
                              ctx,
                              condition,
                              prefix,
                              default_disable_expr);
          }
          break;
        }
      }

      // Get the condition expression(s). More than one `&&&`-joined
      // condition is legal (LRM 12.4.4); AND together the boolean
      // reduction of each condition rather than reading only
      // conditions[0]. A `matches` pattern on any condition introduces
      // destructuring bind semantics this encoder doesn't implement --
      // throw rather than silently evaluate only the plain boolean part
      // of it (mirrors ConditionalExpression handling in
      // expr_encoder.cpp). Each condition's nonzero-reduction is a
      // Bool-sorted term (like LogicalAnd/LogicalOr above use), so the
      // whole conjunction stays Bool.
      Term bool_cond;
      for (auto & c : cond_stmt.conditions) {
        if (c.pattern) {
          throw PonoException(
              "SystemVerilogEncoder: pattern-matching if-condition "
              "('... matches ...') is not supported");
        }
        Term c_bool = expr_encoder_.expr_to_bool(*c.expr, prefix);
        bool_cond =
            bool_cond ? solver_->make_term(And, bool_cond, c_bool) : c_bool;
      }

      // Build then-condition and else-condition.
      Term not_cond = solver_->make_term(Not, bool_cond);
      Term then_cond;
      Term else_cond;
      if (condition == solver_->make_term(true)) {
        // If the outer condition is trivially true, the condition is
        // just the if-expression.
        then_cond = bool_cond;
        else_cond = not_cond;
      } else {
        then_cond = solver_->make_term(And, condition, bool_cond);
        else_cond = solver_->make_term(And, condition, not_cond);
      }

      process_statement(
          cond_stmt.ifTrue, ctx, then_cond, prefix, default_disable_expr);
      if (cond_stmt.ifFalse) {
        process_statement(
            *cond_stmt.ifFalse, ctx, else_cond, prefix, default_disable_expr);
      }
      break;
    }

    case StatementKind::PatternCase: {
      // `case (x) matches` tests each item's *pattern* rather than
      // comparing values, and a pattern can bind names the arm then
      // reads. Unlike the plain Case below, the arms are guarded
      // first-match-wins: a pattern can match anything at all (`.v`
      // does), so without it a later arm would run alongside the one
      // that actually matched.
      auto & pc = stmt.as<PatternCaseStatement>();
      if (pc.condition != CaseStatementCondition::Normal) {
        throw PonoException(
            "SystemVerilogEncoder: only a plain `case ... matches` is "
            "supported; the wildcard and `inside` forms compare pattern "
            "bits in ways a structural match does not");
      }
      Term sel = expr_encoder_.expr_to_term(pc.expr, prefix);
      Term true_term = solver_->make_term(true);
      Term earlier_matched;
      auto & bound = symbol_table_.loop_var_terms();
      for (auto & item : pc.items) {
        std::vector<std::pair<const Symbol *, Term>> bindings;
        Term match = pattern_match(*item.pattern, sel, prefix, bindings);
        if (item.filter) {
          // `matches ... &&& expr`: an extra guard, which may read
          // the names the pattern just bound.
          std::vector<std::pair<const Symbol *, Term>> saved;
          for (auto & b : bindings) {
            auto it = bound.find(b.first);
            saved.emplace_back(b.first,
                               it == bound.end() ? Term() : it->second);
            bound[b.first] = b.second;
          }
          Term guard = expr_encoder_.expr_to_bool(*item.filter, prefix);
          for (auto & sv : saved) {
            if (sv.second) {
              bound[sv.first] = sv.second;
            } else {
              bound.erase(sv.first);
            }
          }
          match = solver_->make_term(And, match, guard);
        }
        Term arm_cond =
            earlier_matched
                ? solver_->make_term(
                      And, match, solver_->make_term(Not, earlier_matched))
                : match;
        earlier_matched = earlier_matched
                              ? solver_->make_term(Or, earlier_matched, match)
                              : match;
        Term full_cond = (condition == true_term)
                             ? arm_cond
                             : solver_->make_term(And, condition, arm_cond);

        std::vector<std::pair<const Symbol *, Term>> saved;
        for (auto & b : bindings) {
          auto it = bound.find(b.first);
          saved.emplace_back(b.first, it == bound.end() ? Term() : it->second);
          bound[b.first] = b.second;
        }
        process_statement(
            *item.stmt, ctx, full_cond, prefix, default_disable_expr);
        for (auto & sv : saved) {
          if (sv.second) {
            bound[sv.first] = sv.second;
          } else {
            bound.erase(sv.first);
          }
        }
      }
      if (pc.defaultCase) {
        Term not_matched = earlier_matched
                               ? solver_->make_term(Not, earlier_matched)
                               : true_term;
        Term default_cond =
            (condition == true_term)
                ? not_matched
                : solver_->make_term(And, condition, not_matched);
        process_statement(
            *pc.defaultCase, ctx, default_cond, prefix, default_disable_expr);
      }
      break;
    }

    case StatementKind::Case: {
      auto & case_stmt = stmt.as<CaseStatement>();
      Term sel = expr_encoder_.expr_to_term(case_stmt.expr, prefix);

      // For casex/casez, a constant item pattern's X (casex only) or
      // Z (both; `?` is just an alias for `z` here) bits are
      // wildcards: build a (mask, value) pair with a 0 at each
      // wildcard bit position and 1s everywhere else, so
      // (sel & mask) == value ignores exactly those positions.
      // Returns nullopt only for a non-constant pattern, where there
      // is nothing to mask against and the caller falls back to plain
      // equality. The pair is MSB-first bit strings rather than
      // integers so that no pattern is too wide to mask: falling back
      // on width alone would hand a pattern with unknown bits to the
      // ordinary literal path, which reads them as unconstrained bits
      // that still have to match rather than as wildcards.
      auto casex_mask = [&](const Expression & pat_expr)
          -> std::optional<std::pair<string, string>> {
        auto cv = pat_expr.eval(expr_encoder_.eval_ctx());
        if (!cv.isInteger()) return std::nullopt;
        auto & sv = cv.integer();
        uint64_t w = pat_expr.type->getBitWidth();
        if (w == 0) return std::nullopt;
        string mask_bits(w, '0'), value_bits(w, '0');
        for (uint64_t i = 0; i < w; ++i) {
          slang::logic_t bit = sv[static_cast<int32_t>(i)];
          bool wildcard =
              (case_stmt.condition == CaseStatementCondition::WildcardXOrZ
               || case_stmt.condition == CaseStatementCondition::Inside)
                  ? bit.isUnknown()
                  : bit.value == slang::logic_t::z.value;
          if (!wildcard) {
            mask_bits[w - 1 - i] = '1';
            if (bit.value == 1) value_bits[w - 1 - i] = '1';
          }
        }
        return std::make_pair(mask_bits, value_bits);
      };
      // `case ... inside` matches by set membership, which for an
      // integral selector is wildcard equality (LRM 11.4.13) -- the
      // same masking casex needs, over both x and z.
      bool is_inside_case =
          case_stmt.condition == CaseStatementCondition::Inside;
      bool is_wildcard_case =
          case_stmt.condition == CaseStatementCondition::WildcardXOrZ
          || case_stmt.condition == CaseStatementCondition::WildcardJustZ
          || is_inside_case;

      Term any_item_matched;
      for (auto & item : case_stmt.items) {
        // Build OR of all patterns matching this item.
        Term item_cond;
        for (auto expr : item.expressions) {
          Term match;
          if (expr->kind == ExpressionKind::ValueRange) {
            // `case (x) inside [lo:hi]:` -- a range rather than a
            // value, which only set membership admits.
            item_cond =
                item_cond ? solver_->make_term(
                                Or,
                                item_cond,
                                expr_encoder_.inside_match(sel, *expr, prefix))
                          : expr_encoder_.inside_match(sel, *expr, prefix);
            continue;
          }
          auto mv = is_wildcard_case ? casex_mask(*expr) : std::nullopt;
          if (mv) {
            uint64_t pat_w = expr->type->getBitWidth();
            Sort pat_sort = solver_->make_sort(BV, pat_w);
            // Raw bit-pattern masks, not numeric values -- always
            // zero-extend.
            Term mask_term =
                resize_to(solver_,
                          solver_->make_term(mv->first, pat_sort, 2),
                          sel->get_sort()->get_width(),
                          false);
            Term value_term =
                resize_to(solver_,
                          solver_->make_term(mv->second, pat_sort, 2),
                          sel->get_sort()->get_width(),
                          false);
            match = solver_->make_term(
                Equal, solver_->make_term(BVAnd, sel, mask_term), value_term);
          } else {
            Term pat = expr_encoder_.expr_to_term(*expr, prefix);
            pat = resize_to(solver_,
                            pat,
                            sel->get_sort()->get_width(),
                            expr->type->isSigned());
            match = solver_->make_term(Equal, sel, pat);
          }
          item_cond =
              item_cond ? solver_->make_term(Or, item_cond, match) : match;
        }
        any_item_matched =
            any_item_matched
                ? solver_->make_term(Or, any_item_matched, item_cond)
                : item_cond;
        Term full_cond = (condition == solver_->make_term(true))
                             ? item_cond
                             : solver_->make_term(And, condition, item_cond);
        process_statement(
            *item.stmt, ctx, full_cond, prefix, default_disable_expr);
      }
      if (case_stmt.defaultCase) {
        // Default: only when none of the other items matched.
        Term not_matched = any_item_matched
                               ? solver_->make_term(Not, any_item_matched)
                               : solver_->make_term(true);
        Term default_cond =
            (condition == solver_->make_term(true))
                ? not_matched
                : solver_->make_term(And, condition, not_matched);
        process_statement(*case_stmt.defaultCase,
                          ctx,
                          default_cond,
                          prefix,
                          default_disable_expr);
      }
      break;
    }

    case StatementKind::Timed: {
      // Skip timing control (e.g., @(posedge clk)) and process the body.
      auto & timed = stmt.as<TimedStatement>();
      process_statement(
          timed.stmt, ctx, condition, prefix, default_disable_expr);
      break;
    }

    case StatementKind::ConcurrentAssertion: {
      auto & ca = stmt.as<ConcurrentAssertionStatement>();
      assertion_walker_.process_concurrent_assertion(
          ca, stmt, prefix, default_disable_expr);
      break;
    }

    case StatementKind::ImmediateAssertion: {
      auto & ia = stmt.as<ImmediateAssertionStatement>();
      assertion_walker_.process_immediate_assertion(ia, condition, prefix);
      break;
    }

    case StatementKind::VariableDeclaration: {
      // Procedural local variable (`int x = ...`): evaluate the
      // initializer once and bind it in the slang EvalContext and our
      // SMT-side loop_var_terms_ map.  A later plain-assignment or
      // `++`/`--` write to it is handled by the ExpressionStatement
      // local-variable fast path above, which re-evaluates the write
      // via slang's constant evaluator and refreshes loop_var_terms_.
      auto & vds = stmt.as<VariableDeclStatement>();
      auto & sym = vds.symbol;
      const Expression * init = sym.getInitializer();
      slang::ConstantValue cv;
      if (init) cv = init->eval(expr_encoder_.eval_ctx());
      if (cv.bad() && !init) cv = sym.getType().getDefaultValue();

      // A constant initial value (or a 2-state type's zero default)
      // also binds the slang-side local, so the unrolling machinery
      // can keep folding this variable in loop bounds and conditions.
      if (cv.isInteger() && !cv.integer().hasUnknown()) {
        expr_encoder_.eval_ctx().createLocal(&sym, cv);
        auto svint = cv.integer();
        uint64_t width = sym.getType().getBitWidth();
        if (width == 0) width = svint.getBitWidth();
        if (width == 0) width = 32;
        Sort sort = solver_->make_sort(BV, width);
        svint.setSigned(false);
        string val_str = svint.toString(slang::LiteralBase::Decimal, false);
        symbol_table_.loop_var_terms()[&sym] =
            solver_->make_term(val_str, sort, 10);
        break;
      }
      // An initializer that isn't elaboration-time constant is still
      // perfectly good logic; bind the term it computes.
      if (init) {
        symbol_table_.loop_var_terms()[&sym] =
            expr_encoder_.expr_to_term(*init, prefix);
        break;
      }
      // No initializer: a 4-state local's default is all-X, which
      // this encoder's 2-valued model cannot represent and must not
      // invent a number for. Leave the variable unbound -- the write
      // that gives it a value binds it, and a read before then is
      // reported by lookup_symbol().
      break;
    }

    case StatementKind::ForLoop: {
      // Compile-time unroll the loop.  The initializers, stop
      // expression, and step expressions are evaluated via
      // expr_encoder_.eval_ctx() below and must succeed as compile-time
      // constants; a non-constant bound throws a PonoException rather than
      // being silently accepted.
      auto & loop = stmt.as<ForLoopStatement>();
      std::vector<const ValueSymbol *> declared;

      auto bind_var = [&](const VariableSymbol & lv) {
        slang::ConstantValue cv;
        if (auto * init = lv.getInitializer()) {
          cv = init->eval(expr_encoder_.eval_ctx());
        }
        if (cv.bad()) {
          cv = lv.getType().getDefaultValue();
        }
        if (!cv.isInteger()) {
          throw PonoException("SystemVerilogEncoder: non-integer for-loop var '"
                              + string(lv.name) + "'");
        }
        expr_encoder_.eval_ctx().createLocal(&lv, cv);
        declared.push_back(&lv);
      };

      auto refresh_bv = [&](const VariableSymbol & lv) {
        auto * cur = expr_encoder_.eval_ctx().findLocal(&lv);
        if (!cur || !cur->isInteger()) {
          throw PonoException("SystemVerilogEncoder: for-loop var '"
                              + string(lv.name) + "' lost its constant value");
        }
        auto svint = cur->integer();
        uint64_t width = lv.getType().getBitWidth();
        if (width == 0) width = svint.getBitWidth();
        if (width == 0) width = 32;
        Sort sort = solver_->make_sort(BV, width);
        svint.setSigned(false);
        string val_str = svint.toString(slang::LiteralBase::Decimal, false);
        symbol_table_.loop_var_terms()[&lv] =
            solver_->make_term(val_str, sort, 10);
      };

      for (auto * lv : loop.loopVars) bind_var(*lv);
      for (auto * init : loop.initializers) {
        if (init->eval(expr_encoder_.eval_ctx()).bad()) {
          throw PonoException(
              "SystemVerilogEncoder: for-loop initializer eval failed");
        }
      }

      constexpr size_t MAX_ITERS = 65536;
      for (size_t it = 0;; ++it) {
        if (it >= MAX_ITERS) {
          throw PonoException("SystemVerilogEncoder: for-loop exceeded "
                              + std::to_string(MAX_ITERS) + " iterations");
        }
        if (loop.stopExpr) {
          auto sv = loop.stopExpr->eval(expr_encoder_.eval_ctx());
          if (sv.bad()) {
            throw PonoException(
                "SystemVerilogEncoder: for-loop stop eval failed");
          }
          if (!sv.isTrue()) break;
        }
        for (auto * lv : loop.loopVars) refresh_bv(*lv);
        bool broke = false;
        try {
          process_statement(
              loop.body, ctx, condition, prefix, default_disable_expr);
        }
        catch (const LoopControlSignal & sig) {
          if (sig.kind == LoopControlSignal::Break) {
            broke = true;
          } else if (sig.kind != LoopControlSignal::Continue) {
            // A Disable targeting some other (typically outer) named
            // block keeps propagating past this loop.
            throw;
          }
          // Continue: swallow and fall through to run the step
          // expressions below, matching SV's `continue` (which still
          // runs the step before the next iteration test) -- same as
          // a normally-completed iteration.
        }
        if (broke) break;
        for (auto * step : loop.steps) {
          if (step->eval(expr_encoder_.eval_ctx()).bad()) {
            throw PonoException(
                "SystemVerilogEncoder: for-loop step eval failed");
          }
        }
      }

      for (auto * sym : declared) {
        symbol_table_.loop_var_terms().erase(sym);
        expr_encoder_.eval_ctx().deleteLocal(sym);
      }
      break;
    }

    case StatementKind::ForeverLoop: {
      // `initial forever @(...) body` (a legacy structural spelling of
      // `always @(...) body`) is recognized and redirected before
      // ever reaching process_statement() at all -- see
      // as_forever_event_body() in process_assignments()/
      // process_instance(). Any `forever` reached *here* is some other
      // shape (no event control, nested inside another statement,
      // etc.): it has no static iteration bound at all, unlike
      // `for`/`while`/`repeat`, which this encoder unrolls up to a
      // compile-time-computable count -- a genuine architectural
      // boundary, not a "not implemented yet" gap. Throw a clear error
      // rather than silently dropping whatever is inside it.
      throw PonoException(
          "SystemVerilogEncoder: 'forever' is only supported as "
          "'initial forever @(...) ...', a structural spelling of "
          "'always @(...) ...'; a bare forever loop has no static "
          "iteration bound and is not supported");
    }

    case StatementKind::WhileLoop: {
      // Compile-time unroll, same contract as `for`: the condition
      // must be constant-evaluable on every iteration (it can
      // reference `for`-loop counters or other locals kept in sync by
      // the ExpressionStatement fast path above) -- a condition that
      // genuinely depends on a runtime (free/registered) signal can't
      // be modeled as C++-level control flow at all, the same
      // architectural boundary as the runtime-dependent break/
      // continue/disable case below.
      auto & loop = stmt.as<WhileLoopStatement>();
      constexpr size_t MAX_ITERS = 65536;
      for (size_t it = 0;; ++it) {
        auto cv = loop.cond.eval(expr_encoder_.eval_ctx());
        if (cv.bad()) {
          throw PonoException(
              "SystemVerilogEncoder: 'while' condition is not a "
              "compile-time constant (runtime-dependent while loops "
              "are not supported)");
        }
        if (!cv.isTrue()) break;
        if (it >= MAX_ITERS) {
          throw PonoException("SystemVerilogEncoder: 'while' loop exceeded "
                              + std::to_string(MAX_ITERS) + " iterations");
        }
        bool broke = false;
        try {
          process_statement(
              loop.body, ctx, condition, prefix, default_disable_expr);
        }
        catch (const LoopControlSignal & sig) {
          if (sig.kind == LoopControlSignal::Break) {
            broke = true;
          } else if (sig.kind != LoopControlSignal::Continue) {
            throw;
          }
        }
        if (broke) break;
      }
      break;
    }

    case StatementKind::DoWhileLoop: {
      // Same as WhileLoop, but the condition is tested after the
      // first execution of the body (`do ... while (cond);`).
      auto & loop = stmt.as<DoWhileLoopStatement>();
      constexpr size_t MAX_ITERS = 65536;
      for (size_t it = 0;; ++it) {
        if (it >= MAX_ITERS) {
          throw PonoException("SystemVerilogEncoder: 'do-while' loop exceeded "
                              + std::to_string(MAX_ITERS) + " iterations");
        }
        bool broke = false;
        try {
          process_statement(
              loop.body, ctx, condition, prefix, default_disable_expr);
        }
        catch (const LoopControlSignal & sig) {
          if (sig.kind == LoopControlSignal::Break) {
            broke = true;
          } else if (sig.kind != LoopControlSignal::Continue) {
            throw;
          }
        }
        if (broke) break;
        auto cv = loop.cond.eval(expr_encoder_.eval_ctx());
        if (cv.bad()) {
          throw PonoException(
              "SystemVerilogEncoder: 'do-while' condition is not a "
              "compile-time constant (runtime-dependent do-while "
              "loops are not supported)");
        }
        if (!cv.isTrue()) break;
      }
      break;
    }

    case StatementKind::RepeatLoop: {
      // The trip count is evaluated once, up front (per the LRM,
      // `repeat` takes a plain expression, not a re-tested condition
      // like `while`); an unroll-time-unresolvable count (a runtime
      // signal) is out of scope, same contract as `for`/`while`
      // bounds.
      auto & loop = stmt.as<RepeatLoopStatement>();
      auto cv = loop.count.eval(expr_encoder_.eval_ctx());
      auto n_opt = cv.bad() ? std::nullopt : cv.integer().as<uint64_t>();
      if (!n_opt) {
        throw PonoException(
            "SystemVerilogEncoder: 'repeat' count is not a "
            "compile-time constant (runtime-dependent repeat counts "
            "are not supported)");
      }
      constexpr uint64_t MAX_ITERS = 65536;
      uint64_t n = *n_opt;
      if (n > MAX_ITERS) {
        throw PonoException("SystemVerilogEncoder: 'repeat' loop exceeded "
                            + std::to_string(MAX_ITERS) + " iterations");
      }
      for (uint64_t it = 0; it < n; ++it) {
        bool broke = false;
        try {
          process_statement(
              loop.body, ctx, condition, prefix, default_disable_expr);
        }
        catch (const LoopControlSignal & sig) {
          if (sig.kind == LoopControlSignal::Break) {
            broke = true;
          } else if (sig.kind != LoopControlSignal::Continue) {
            throw;
          }
        }
        if (broke) break;
      }
      break;
    }

    case StatementKind::ForeachLoop: {
      // Scoped to a single iterated dimension with a concrete loop
      // variable and a statically-known range -- exactly the shape
      // `foreach (arr[i])` produces for a fixed-size packed
      // array/vector, which is all a compile-time-unrolling model can
      // support. Multiple dimensions (`foreach (arr[i][j])`) or a
      // dynamically-sized dimension (a real dynamic array/queue,
      // which this encoder has no sort for anyway) throw a clear
      // error instead of silently iterating the wrong thing or just
      // the first dimension.
      auto & loop = stmt.as<ForeachLoopStatement>();
      if (loop.loopDims.size() != 1 || !loop.loopDims[0].loopVar
          || !loop.loopDims[0].range) {
        throw PonoException(
            "SystemVerilogEncoder: 'foreach' is only supported over a "
            "single statically-sized dimension with a loop variable "
            "(multi-dimensional or dynamically-sized foreach is not "
            "supported)");
      }
      auto & dim = loop.loopDims[0];
      auto & iter_sym = *dim.loopVar;
      int32_t lo = dim.range->lower();
      int32_t hi = dim.range->upper();
      uint64_t width = iter_sym.getType().getBitWidth();
      if (width == 0) width = 32;

      constexpr uint64_t MAX_ITERS = 65536;
      if (static_cast<uint64_t>(hi) - static_cast<uint64_t>(lo) + 1
          > MAX_ITERS) {
        throw PonoException("SystemVerilogEncoder: 'foreach' loop exceeded "
                            + std::to_string(MAX_ITERS) + " iterations");
      }
      for (int32_t idx = lo; idx <= hi; ++idx) {
        slang::SVInt iv(static_cast<slang::bitwidth_t>(width),
                        static_cast<uint64_t>(idx),
                        /*isSigned=*/true);
        expr_encoder_.eval_ctx().createLocal(&iter_sym,
                                             slang::ConstantValue(iv));
        refresh_loop_var_term(iter_sym);
        bool broke = false;
        try {
          process_statement(
              loop.body, ctx, condition, prefix, default_disable_expr);
        }
        catch (const LoopControlSignal & sig) {
          if (sig.kind != LoopControlSignal::Break
              && sig.kind != LoopControlSignal::Continue) {
            symbol_table_.loop_var_terms().erase(&iter_sym);
            expr_encoder_.eval_ctx().deleteLocal(&iter_sym);
            throw;
          }
          broke = sig.kind == LoopControlSignal::Break;
        }
        symbol_table_.loop_var_terms().erase(&iter_sym);
        expr_encoder_.eval_ctx().deleteLocal(&iter_sym);
        if (broke) break;
      }
      break;
    }

    case StatementKind::Return: {
      auto & ret = stmt.as<ReturnStatement>();
      if (!current_return_var_) {
        throw PonoException(
            "SystemVerilogEncoder: `return` outside an inlined subroutine");
      }
      // A void `return;` ends the call without producing a value.
      if (!ret.expr) break;
      if (condition != solver_->make_term(true)) {
        // Statements after this one would still be walked, so the
        // value would be whatever the last reached assignment left,
        // not this one. Assigning the function name under the
        // condition instead expresses the same thing and is modelled.
        throw PonoException(
            "SystemVerilogEncoder: a `return` reached only under a "
            "runtime condition is not supported; assign the function's "
            "name instead");
      }
      symbol_table_.loop_var_terms()[current_return_var_] =
          expr_encoder_.expr_to_term(*ret.expr, prefix);
      break;
    }

    case StatementKind::Break:
    case StatementKind::Continue:
    case StatementKind::Disable: {
      // `condition` is only ever narrowed away from the trivial `true`
      // term by the Conditional/Case cases' *general symbolic* guard
      // building -- the Conditional case's constant-fold fast path
      // (see its comment) always passes `condition` through unchanged,
      // and Block/ForLoop never touch it at all. So reaching this
      // statement with `condition` still exactly `true` means every
      // enclosing `if` along the way was compile-time-resolved, i.e.
      // this control-flow statement is genuinely, unconditionally
      // reached at this point in the unrolling -- interpretable via
      // LoopControlSignal, caught by the nearest enclosing loop
      // (Break/Continue) or matching named Block (Disable). Any other
      // value means we came through at least one runtime-dependent
      // `if`, which (since that path processes both arms
      // unconditionally) can't be correctly modeled as C++-level
      // control flow at all -- a clear error beats silently always- or
      // never-triggering regardless of the real condition.
      if (condition != solver_->make_term(true)) {
        throw PonoException(
            "SystemVerilogEncoder: break/continue/disable is only "
            "supported when its controlling condition is a compile-time "
            "constant (e.g. depends only on already-unrolled for-loop "
            "counters)");
      }
      if (stmt.kind == StatementKind::Break) {
        throw LoopControlSignal{ LoopControlSignal::Break };
      }
      if (stmt.kind == StatementKind::Continue) {
        throw LoopControlSignal{ LoopControlSignal::Continue };
      }
      auto & ds = stmt.as<DisableStatement>();
      const Symbol * target = nullptr;
      if (ds.target.kind == ExpressionKind::ArbitrarySymbol) {
        target = ds.target.as<ArbitrarySymbolExpression>().symbol.get();
      }
      throw LoopControlSignal{ LoopControlSignal::Disable, target };
    }

    default:
      // Other statement kinds (e.g. wait, return, event triggers,
      // randcase): not supported in synthesizable subset. Log a
      // warning and skip.
      logger.log(1,
                 "SystemVerilogEncoder: skipping unsupported statement kind {}",
                 static_cast<int>(stmt.kind));
      break;
  }
}

}  // namespace pono
