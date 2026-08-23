/*! \file ast_helpers.h
 *  \brief Shared free-function AST helpers used across multiple
 *         SystemVerilogEncoder translation units.
 */
#pragma once

#include <cstdint>
#include <optional>
#include <unordered_set>

#include "slang/ast/EvalContext.h"
#include "slang/ast/Expression.h"
#include "slang/ast/Statement.h"
#include "slang/ast/Symbol.h"

namespace pono {

// Helper to iterate over sub-statements of a BlockStatement body.
// The body is a single Statement; if it is a StatementList, iterate its
// children, otherwise visit the single statement directly.
template <typename Func>
void for_each_stmt_in_block(const slang::ast::BlockStatement & block,
                            Func && func)
{
  auto & body = block.body;
  if (body.kind == slang::ast::StatementKind::List) {
    auto & list = body.as<slang::ast::StatementList>();
    for (auto * s : list.list) {
      func(*s);
    }
  } else {
    func(body);
  }
}

// Collects targets of non-blocking assignments inside `stmt` (recursing
// through blocks/conditionals/case/loops). Used both by pre-scan (to
// classify state vars) and by process_assignments()/process_instance()
// (to decide whether a wire's driving process uses blocking or
// non-blocking writes).
void collect_nonblocking_targets(
    const slang::ast::Statement & stmt,
    std::unordered_set<const slang::ast::Symbol *> & targets);

// A modport-qualified interface port access (`b.data`, where `b`'s
// declared type is e.g. `bus_if.master`) resolves to a synthesized
// ModportPortSymbol proxy, not directly to the interface instance's
// own `data` VariableSymbol -- even though the plain (non-modport)
// interface case (`b.data` where `b`'s type is just `bus_if`) resolves
// straight to that same shared symbol. Redirect through
// ModportPortSymbol::internalSymbol so every access path -- with or
// without a modport qualifier -- converges on the same underlying
// symbol identity. A no-op for every other symbol kind.
const slang::ast::Symbol & canonicalize_modport_port(
    const slang::ast::Symbol & sym);

// Identifies the base ValueSymbol underlying a (possibly nested)
// bit/range-select or struct-member-access LHS.  Returns nullptr if
// the LHS shape isn't supported by the encoder (e.g. concatenation LHS).
const slang::ast::Symbol * find_lhs_base(const slang::ast::Expression & lhs);

// Describes an LHS slice: which base symbol gets written and at what
// bit range.  For a NamedValue (or HierarchicalValue) LHS this is the
// full range [0, base_w-1]; for nested ElementSelects/RangeSelects of
// constant indices/bounds, or MemberAccess field selects, the range
// narrows accordingly while base_w stays the full base bit width.
struct LValueDesc
{
  const slang::ast::Symbol * base;
  uint64_t lo;
  uint64_t hi;
  uint64_t base_w;
};

// resolve_lvalue() has exactly one legitimate reason to return nullopt
// rather than throw: ExpressionKind::ElementSelect with a genuinely
// non-constant (runtime-variable) index, which the caller detects via
// `lhs_expr.kind == ExpressionKind::ElementSelect` and re-dispatches to
// process_dynamic_element_assign() instead. Every other unresolvable
// shape -- a whole ExpressionKind this function has no case for at
// all, or a malformed/out-of-bounds constant within a case it does
// recognize -- throws instead of silently dropping the write, per the
// same "throw rather than silently mis-encode" contract enforced
// everywhere else in this file (see expr_to_term()'s own default
// case). A nested dynamic index (e.g. `arr[i][3:0] <= ...`) still
// propagates as nullopt through `inner`, since the top-level caller's
// ElementSelect-shaped fallback only ever re-dispatches on the
// outermost expression.
std::optional<LValueDesc> resolve_lvalue(const slang::ast::Expression & lhs,
                                         slang::ast::EvalContext & ctx);

// Internal control-flow signal for `break`/`continue`/`disable`,
// thrown by process_statement()'s Break/Continue/Disable cases and
// caught by whichever enclosing construct can absorb it: any of the
// unrolled loop kinds -- ForLoop, WhileLoop, DoWhileLoop, RepeatLoop,
// ForeachLoop -- (for Break/Continue) or a named Block whose symbol
// matches disable_target (for Disable). This only correctly models
// compile-time-reachable control flow -- e.g. `if (i == 2) break;`
// where `i` is an already-unrolled `for`-loop counter, handled by the
// Conditional case's constant-fold fast path, which branches in C++
// rather than building a symbolic guard for both arms, so this signal
// can propagate out of only the taken branch exactly like a real
// break/continue/disable would.  A signal thrown from inside a
// *runtime*-dependent condition (which the general symbolic-guard
// path processes unconditionally for both arms) has no correct
// interpretation here and is never meant to be caught by anything;
// process_always_ff()/process_always_comb()/process_initial() convert
// any instance that escapes all the way out into a clear
// PonoException rather than letting it propagate as a raw internal
// exception type.
struct LoopControlSignal
{
  enum Kind
  {
    Break,
    Continue,
    Disable
  } kind;
  const slang::ast::Symbol * disable_target = nullptr;  // only for Disable
};

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
    const slang::ast::Statement & stmt);

}  // namespace pono
