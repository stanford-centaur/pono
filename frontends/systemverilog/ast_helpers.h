/*!
 * \file ast_helpers.h
 * \brief LHS-resolution, traversal, and control-flow helper declarations.
 * \author Áron Ricardo Perez-Lopez
 * \date 2026
 * \copyright See the LICENSE file in the top-level source directory.
 *
 * These helpers fall into three groups. Statement traversal
 * (for_each_stmt_in_block, collect_nonblocking_targets) walks
 * block/conditional/case/loop bodies to find non-blocking-assignment targets
 * during pre-scan and process_instance(). LHS/lvalue resolution
 * (canonicalize_signal_alias, find_lhs_base, resolve_lvalue, LValueDesc) maps
 * an assignment's left-hand side through modport indirection down to a base
 * Symbol and constant bit range. Compile-time control flow (LoopControlSignal,
 * as_forever_event_body) lets process_statement() model break/continue/disable
 * across unrolled loops and recognize the legacy `initial forever @(...) body`
 * idiom as equivalent to an always block.
 */
#pragma once

#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <unordered_set>

#include "slang/ast/EvalContext.h"
#include "slang/ast/Expression.h"
#include "slang/ast/Statement.h"
#include "slang/ast/Symbol.h"

namespace slang::ast {
class Scope;
class ElementSelectExpression;
class PortSymbol;
}  // namespace slang::ast

namespace pono {

/** Invoke `fn` for every concrete member of `scope`, recursing through
 *  generate-block scopes so the caller sees the unrolled per-iteration
 *  members directly. Uninstantiated generate blocks (the unselected arm
 *  of a generate-if / generate-case) are skipped.
 *
 *  Takes a type-erased `std::function` (rather than a template
 *  parameter) specifically so its single definition can live in exactly
 *  one file while still being callable, as an ordinary non-template
 *  function, from every other file in this directory -- a template's
 *  definition would instead need to be visible (and separately
 *  instantiated) in each of those translation units, which for this
 *  function would also force this header to fully include (rather than
 *  just forward-declare) the slang types its body depends on.
 *  @param scope the scope whose members should be visited
 *  @param prefix the caller's current hierarchical name prefix --
 *         updated in place while descending into a generate-for/
 *         instance-array's bracket-indexed child scopes, and restored
 *         before returning
 *  @param fn callback invoked once per concrete member symbol
 */
void walk_members(const slang::ast::Scope & scope,
                  std::string & prefix,
                  const std::function<void(const slang::ast::Symbol &)> & fn);

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
const slang::ast::Symbol & canonicalize_signal_alias(
    const slang::ast::Symbol & sym);

// The instance-internal symbol a port connects to, which is what
// every part of the encoder keys a port on: the declaration that
// gives it a term, and the instance connection that binds it to the
// parent side. Usually PortSymbol::internalSymbol, which slang fills
// in for every ordinary port. An explicit port (`output .o(w)`)
// names its internal signal through an expression instead and leaves
// that field null; so does an empty slot in a port list, which
// connects to nothing at all and returns nullptr here.
// Throws for an explicit port bound to anything but a whole signal:
// one bound to a select or a concatenation stands for part of a
// symbol, or parts of several, and there is no single term for the
// two sides to share.
const slang::ast::Symbol * port_internal_symbol(
    const slang::ast::PortSymbol & port);

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
//
// An element of an *unpacked* array is not a bit range of its base, so
// LValueDesc cannot name it and a caller that passes no `array_elem`
// gets nullopt. Passing one instead makes that element a synthetic
// base: `*array_elem` is set to the select naming it, `base` comes
// back nullptr, and `lo`/`hi`/`base_w` describe the written range
// *within the element*, so the MemberAccess/RangeSelect/ElementSelect
// layers above it (`mem[i].f`, `mem[i][3:0]`, `mem[i][j]`) compose
// their offsets exactly as they do over a real symbol.
std::optional<LValueDesc> resolve_lvalue(
    const slang::ast::Expression & lhs,
    slang::ast::EvalContext & ctx,
    const slang::ast::ElementSelectExpression ** array_elem = nullptr);

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

// True if `body` is guarded by an edge-sensitive event control
// (`@(posedge clk)`, `@(negedge rst_n)`, or an event list containing
// one), i.e. the block is a register rather than combinational logic.
// `always_ff` says so in its keyword, but a plain `always` does not --
// only its timing control distinguishes `always @(posedge clk)` from
// `always @(*)`, and the two must reach opposite halves of
// process_instance()'s walk.
bool is_edge_triggered(const slang::ast::Statement & body);

/** Whether `body` consists of concurrent assertions and nothing else
 *  (an empty block counts). Used to pick the parts of a `program`
 *  worth encoding: its stimulus is simulation-only, but an assertion
 *  written inside it is an ordinary property.
 */
bool is_concurrent_assertion_only(const slang::ast::Statement & body);

// Collects the block-locals of `body` that some execution path reads
// before writing. Those are not temporaries at all: a variable read
// where no path assigned it keeps its previous value, which is
// storage (a flop in a clocked block, a latch in a combinational
// one). Every other local is bound to the term its write computes
// and needs no state.
//
// Conservative in the safe direction -- a branch contributes only
// what all of its arms assign, and a loop body contributes nothing
// since it may run zero times -- so it can name a local that does
// not really need storage, never miss one that does.
/** The symbols `body` assigns on *every* path through it. A branch
 *  contributes only what all of its arms assign, and a loop body
 *  nothing, since it may run zero times. Anything a block assigns
 *  but that is missing here keeps its old value on some path --
 *  which in a combinational block is a latch.
 */
void collect_definitely_assigned(
    const slang::ast::Statement & body,
    std::unordered_set<const slang::ast::Symbol *> & out);

void collect_hold_locals(const slang::ast::Statement & body,
                         std::unordered_set<const slang::ast::Symbol *> & out);

// True if `sym` belongs to a procedural block or a subroutine rather
// than to module scope: a temporary with no life beyond one execution
// or one call, so it must never be classified as a register.
bool is_block_local(const slang::ast::Symbol & sym);

}  // namespace pono
