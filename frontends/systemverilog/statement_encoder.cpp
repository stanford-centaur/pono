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
#include <string>
#include <vector>

#include "frontends/systemverilog/assertion_walker.h"
#include "frontends/systemverilog/ast_helpers.h"
#include "frontends/systemverilog/bit_utils.h"
#include "frontends/systemverilog/expr_encoder.h"
#include "frontends/systemverilog/symbol_table.h"
#include "slang/ast/EvalContext.h"
#include "slang/ast/Expression.h"
#include "slang/ast/SemanticFacts.h"
#include "slang/ast/Statement.h"
#include "slang/ast/Symbol.h"
#include "slang/ast/expressions/AssignmentExpressions.h"
#include "slang/ast/expressions/MiscExpressions.h"
#include "slang/ast/expressions/OperatorExpressions.h"
#include "slang/ast/expressions/SelectExpressions.h"
#include "slang/ast/statements/ConditionalStatements.h"
#include "slang/ast/statements/LoopStatements.h"
#include "slang/ast/statements/MiscStatements.h"
#include "slang/ast/symbols/BlockSymbols.h"
#include "slang/ast/symbols/VariableSymbols.h"
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

void StatementEncoder::process_dynamic_element_assign(
    const slang::ast::ElementSelectExpression & sel,
    const slang::ast::Expression & rhs_expr,
    StmtContext ctx,
    const Term & condition,
    const string & prefix)
{
  using namespace slang::ast;

  // Only a direct `base[idx] = rhs` pattern is supported: the select
  // must sit directly on a plain variable, not on a nested select.
  if (sel.value().kind != ExpressionKind::NamedValue
      && sel.value().kind != ExpressionKind::HierarchicalValue) {
    return;
  }
  const Symbol * sym = &canonicalize_modport_port(
      (sel.value().kind == ExpressionKind::NamedValue)
          ? sel.value().as<NamedValueExpression>().symbol
          : sel.value().as<HierarchicalValueExpression>().symbol);
  bool aliased = symbol_table_.port_output_aliases().count(sym) > 0;
  {
    uint64_t sym_w = sym->as<ValueSymbol>().getType().getBitWidth();
    auto pieces = symbol_table_.resolve_output_alias_pieces(sym, 0, sym_w - 1);
    uint64_t piece_w =
        pieces.empty()
            ? 0
            : pieces[0].sym->as<ValueSymbol>().getType().getBitWidth();
    if (pieces.size() != 1 || pieces[0].target_lo != 0
        || pieces[0].target_hi + 1 != piece_w) {
      // Dynamic-index writes into a bus-element alias (e.g. one
      // element of an instance array wired to a slice of a parent
      // bus) or a concatenation-target alias aren't supported; leave
      // unconstrained rather than risk a wrong encoding.
      return;
    }
    sym = pieces[0].sym;
  }

  uint64_t elem_w = sel.type->getBitWidth();
  if (elem_w == 0) elem_w = 1;

  bool wire_comb = ctx == StmtContext::COMBINATIONAL
                   && symbol_table_.wire_symbols().count(sym);
  Term prev_base;
  Term state_term;
  if (wire_comb) {
    auto pit = symbol_table_.pending_comb_updates().find(sym);
    if (pit != symbol_table_.pending_comb_updates().end()) {
      prev_base = pit->second;
    } else {
      auto sit = symbol_table_.symbol_to_term().find(sym);
      if (sit != symbol_table_.symbol_to_term().end()) prev_base = sit->second;
    }
  } else if (ctx == StmtContext::NEXT_STATE) {
    auto sit = symbol_table_.symbol_to_term().find(sym);
    if (sit == symbol_table_.symbol_to_term().end()) return;
    state_term = sit->second;
    auto pit = symbol_table_.pending_next_updates().find(state_term);
    prev_base = (pit != symbol_table_.pending_next_updates().end())
                    ? pit->second
                    : state_term;
  } else {
    // COMBINATIONAL non-wire / INITIAL dynamic-index writes aren't
    // needed by any currently-supported construct; leave unsupported
    // rather than risk an under-constrained partial encoding.
    return;
  }
  if (!prev_base) return;

  Term idx = expr_encoder_.expr_to_term(sel.selector(), prefix);
  Term rhs = expr_encoder_.expr_to_term(rhs_expr, prefix);
  rhs = resize_to(solver_, rhs, elem_w, sel.type->isSigned());
  Term combined = replace_bits_dynamic(solver_, prev_base, rhs, idx, elem_w);

  if (wire_comb) {
    if (aliased) symbol_table_.pending_comb_aliased().insert(sym);
    symbol_table_.pending_comb_updates()[sym] =
        (condition == solver_->make_term(true))
            ? combined
            : solver_->make_term(Ite, condition, combined, prev_base);
    return;
  }

  // ctx == NEXT_STATE
  symbol_table_.pending_next_updates()[state_term] =
      (condition == solver_->make_term(true))
          ? combined
          : solver_->make_term(Ite, condition, combined, prev_base);
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
            if (expr_encoder_.eval_ctx().findLocal(&vsym)) {
              if (expr.eval(expr_encoder_.eval_ctx()).bad()) {
                throw PonoException(
                    "SystemVerilogEncoder: assignment to local variable '"
                    + string(base->name) + "' failed to constant-evaluate");
              }
              refresh_loop_var_term(vsym);
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
        if (lhs_expr.kind == ExpressionKind::Concatenation) {
          uint64_t total_w = lhs_expr.type->getBitWidth();
          if (total_w == 0) return {};
          std::vector<LValueWrite> writes;
          uint64_t covered = 0;
          for (auto * operand :
               lhs_expr.as<ConcatenationExpression>().operands()) {
            uint64_t seg_w = operand->type->getBitWidth();
            if (seg_w == 0 || covered + seg_w > total_w) return {};
            auto sub = begin_write(*operand);
            if (sub.empty()) return {};
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
          uint64_t sym_w = sym->as<ValueSymbol>().getType().getBitWidth();
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
            uint64_t sym_w = w.sym->as<ValueSymbol>().getType().getBitWidth();
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
            // Non-wire LHS (e.g. output port reg, or a partially-
            // written base): constrain the appropriate slice via
            // add_constraint, which accepts terms involving input
            // vars (so RHSes that reference input ports work).
            Term lhs_slice =
                full_write
                    ? lhs_term
                    : solver_->make_term(Op(Extract, w.hi, w.lo), lhs_term);
            Term eq = solver_->make_term(Equal, lhs_slice, rhs);
            if (condition != solver_->make_term(true)) {
              eq = solver_->make_term(Implies, condition, eq);
            }
            fts_.add_constraint(eq);
            break;
          }
          case StmtContext::INITIAL: {
            Term lhs_slice =
                full_write
                    ? lhs_term
                    : solver_->make_term(Op(Extract, w.hi, w.lo), lhs_term);
            Term eq = solver_->make_term(Equal, lhs_slice, rhs);
            fts_.constrain_init(eq);
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

        auto writes = begin_write(lhs_expr);
        if (writes.empty()) {
          // resolve_lvalue() only handles constant-index selects; a
          // runtime-variable index (`arr[idx] = rhs`) needs a
          // dynamic-position splice instead of a static bit range.
          if (lhs_expr.kind == ExpressionKind::ElementSelect) {
            process_dynamic_element_assign(
                lhs_expr.as<ElementSelectExpression>(),
                rhs_expr,
                ctx,
                condition,
                prefix);
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
        // A concatenation-target LHS is always unsigned per the LRM
        // (positional bit-splicing, not a numeric value); otherwise
        // use the RHS expression's own signedness.
        bool rhs_signed = lhs_expr.kind != ExpressionKind::Concatenation
                          && rhs_expr.type->isSigned();
        rhs = resize_to(solver_, rhs, total_w, rhs_signed);
        for (auto & w : writes) {
          commit_write(w, slice_of(rhs, w.rhs_lo, w.rhs_hi));
        }
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
      // Bool-sorted term (like LogicalAnd/LogicalOr above use); only
      // convert back to BV1 once, after ANDing, to avoid mixing Bool
      // and BV1 sorts under a single BV-typed `And`.
      Term bool_cond;
      for (auto & c : cond_stmt.conditions) {
        if (c.pattern) {
          throw PonoException(
              "SystemVerilogEncoder: pattern-matching if-condition "
              "('... matches ...') is not supported");
        }
        Term c_term = expr_encoder_.expr_to_term(*c.expr, prefix);
        Term c_bool = solver_->make_term(
            Distinct, c_term, solver_->make_term(0, c_term->get_sort()));
        bool_cond =
            bool_cond ? solver_->make_term(And, bool_cond, c_bool) : c_bool;
      }

      // Build then-condition and else-condition.
      Sort bv1 = solver_->make_sort(BV, 1);
      Term one = solver_->make_term(1, bv1);
      Term zero = solver_->make_term(0, bv1);
      Term cond_term = solver_->make_term(Ite, bool_cond, one, zero);

      Term then_cond;
      Term else_cond;
      if (condition == solver_->make_term(true)) {
        // If the outer condition is trivially true, the condition is
        // just the if-expression.
        then_cond = solver_->make_term(Equal, cond_term, one);
        else_cond = solver_->make_term(Equal, cond_term, zero);
      } else {
        Term cond_eq_one = solver_->make_term(Equal, cond_term, one);
        Term cond_eq_zero = solver_->make_term(Equal, cond_term, zero);
        then_cond = solver_->make_term(And, condition, cond_eq_one);
        else_cond = solver_->make_term(And, condition, cond_eq_zero);
      }

      process_statement(
          cond_stmt.ifTrue, ctx, then_cond, prefix, default_disable_expr);
      if (cond_stmt.ifFalse) {
        process_statement(
            *cond_stmt.ifFalse, ctx, else_cond, prefix, default_disable_expr);
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
      // Returns nullopt for a non-constant pattern (nothing to mask
      // against) or a width too large for this uint64_t-based mask,
      // in which case the caller falls back to plain equality.
      auto casex_mask = [&](const Expression & pat_expr)
          -> std::optional<std::pair<uint64_t, uint64_t>> {
        auto cv = pat_expr.eval(expr_encoder_.eval_ctx());
        if (!cv.isInteger()) return std::nullopt;
        auto & sv = cv.integer();
        uint64_t w = pat_expr.type->getBitWidth();
        if (w == 0 || w > 64) return std::nullopt;
        uint64_t mask = 0, value = 0;
        for (uint64_t i = 0; i < w; ++i) {
          slang::logic_t bit = sv[static_cast<int32_t>(i)];
          bool wildcard =
              (case_stmt.condition == CaseStatementCondition::WildcardXOrZ)
                  ? bit.isUnknown()
                  : bit.value == slang::logic_t::z.value;
          if (!wildcard) {
            mask |= (uint64_t{ 1 } << i);
            if (bit.value == 1) value |= (uint64_t{ 1 } << i);
          }
        }
        return std::make_pair(mask, value);
      };
      bool is_wildcard_case =
          case_stmt.condition == CaseStatementCondition::WildcardXOrZ
          || case_stmt.condition == CaseStatementCondition::WildcardJustZ;

      Term any_item_matched;
      for (auto & item : case_stmt.items) {
        // Build OR of all patterns matching this item.
        Term item_cond;
        for (auto expr : item.expressions) {
          Term match;
          auto mv = is_wildcard_case ? casex_mask(*expr) : std::nullopt;
          if (mv) {
            uint64_t pat_w = expr->type->getBitWidth();
            Sort pat_sort = solver_->make_sort(BV, pat_w);
            // Raw bit-pattern masks, not numeric values -- always
            // zero-extend.
            Term mask_term = resize_to(
                solver_,
                solver_->make_term(std::to_string(mv->first), pat_sort, 10),
                sel->get_sort()->get_width(),
                false);
            Term value_term = resize_to(
                solver_,
                solver_->make_term(std::to_string(mv->second), pat_sort, 10),
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
      slang::ConstantValue cv;
      if (auto * init = sym.getInitializer()) {
        cv = init->eval(expr_encoder_.eval_ctx());
      }
      if (cv.bad()) {
        cv = sym.getType().getDefaultValue();
      }
      if (!cv.isInteger()) {
        throw PonoException("SystemVerilogEncoder: non-integer local '"
                            + string(sym.name) + "'");
      }
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
