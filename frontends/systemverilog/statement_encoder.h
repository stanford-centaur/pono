/*!
 * \file statement_encoder.h
 * \brief The process_statement() switch encoding SV procedural statements.
 * \author Áron Ricardo Perez-Lopez
 * \date 2026
 * \copyright See the LICENSE file in the top-level source directory.
 *
 * StatementEncoder dispatches on slang::ast::StatementKind to encode
 * assignments (plain, compound, ++/--, concatenation-target and dynamic-
 * index LHS splicing), if/case/casex/casez, loops, and concurrent/
 * immediate assertions into the FunctionalTransitionSystem. Loops (for/
 * while/do-while/repeat/foreach) are unrolled at compile time via slang's
 * own constant evaluator, up to a fixed iteration cap; a genuinely
 * runtime-dependent loop bound, or a break/continue/disable reached
 * through a non-constant condition, is rejected rather than silently
 * mis-encoded. Local (non-state, non-wire) variables are mirrored as SMT
 * constants and kept in sync with slang's evaluator by
 * refresh_loop_var_term() after each constant-evaluated write.
 * break/continue/disable are modeled as C++ exceptions (LoopControlSignal)
 * caught by the nearest loop or matching named block. Assertions dispatch
 * to AssertionWalker.
 *
 * Depends on SymbolTable, ExprEncoder, and AssertionWalker -- all already
 * independent of the rest of the encoder, so this class is too: it holds
 * no reference back to SystemVerilogEncoder. The hierarchical name prefix
 * and the enclosing module's `default disable iff` expression (needed
 * only for a ConcurrentAssertion statement, and otherwise threaded
 * through unused) are explicit parameters to process_statement() rather
 * than ambient state, matching SymbolTable's/AssertionWalker's own
 * design; the caller resolves `default disable iff` once per instance
 * body via Compilation::getDefaultDisable(), so this class needs no
 * Compilation/Scope dependency of its own.
 */
#pragma once

#include <string>
#include <unordered_set>

#include "core/fts.h"
#include "frontends/systemverilog/expr_encoder.h"
#include "smt-switch/smt.h"

namespace slang {
class ConstantValue;
}

namespace slang::ast {
class ElementSelectExpression;
class FixedSizeUnpackedArrayType;
class Expression;
class Statement;
class Symbol;
class Type;
class ValueSymbol;
}  // namespace slang::ast

namespace pono {

class AssertionWalker;
class ExprEncoder;
class SymbolTable;

class StatementEncoder : public ExprEncoder::SubroutineInliner
{
 public:
  /** Context for statement processing: whether we are building next-state
   *  updates (always_ff), combinational definitions (always_comb), or
   *  initial constraints.
   */
  enum class StmtContext
  {
    NEXT_STATE,     ///< Inside always_ff (also always_latch and a legacy
                    ///< `forever @(...)` spelling of always_ff): build
                    ///< next-state functions
    COMBINATIONAL,  ///< Inside always_comb: build combinational definitions
    INITIAL         ///< Inside initial: build init constraints
  };

  StatementEncoder(SymbolTable & symbol_table,
                   ExprEncoder & expr_encoder,
                   AssertionWalker & assertion_walker,
                   FunctionalTransitionSystem & fts,
                   const smt::SmtSolver & solver);

  /** Recursively process a statement, extracting assignments.
   *  @param stmt the slang statement to process
   *  @param ctx  what kind of block we are in
   *  @param condition accumulated path condition (for if/case nesting)
   *  @param prefix the caller's current hierarchical name prefix
   *  @param default_disable_expr the enclosing module's `default
   *         disable iff` condition, used only if a ConcurrentAssertion
   *         statement is reached and has no explicit `disable iff` of
   *         its own; null if none applies
   */
  void process_statement(const slang::ast::Statement & stmt,
                         StmtContext ctx,
                         const smt::Term & condition,
                         const std::string & prefix,
                         const slang::ast::Expression * default_disable_expr);

  /** Walk an inlined subroutine body, so a call in an expression can
   *  be given a value. Formals and the return variable are bound by
   *  ExprEncoder::inline_call() before this runs; the body's writes
   *  land on them through the ordinary local-write path.
   */
  void inline_subroutine_body(const slang::ast::Statement & body,
                              const slang::ast::Symbol & return_var,
                              const std::string & prefix) override;

 private:
  /** The same walk for a body that produces no value, which is what
   *  a task call as a statement needs. */
  void inline_subroutine_body_no_return(const slang::ast::Statement & body,
                                        const std::string & prefix);

  /** Handle `base[idx] = rhs` (nonblocking or blocking) when `idx` is
   *  not a compile-time constant, so resolve_lvalue() can't produce a
   *  static bit range.  Only a direct select on a plain variable base
   *  is supported (no nested dynamic selects).  A no-op if the base
   *  isn't a plain variable, if it resolves through
   *  port_output_aliases_ to anything other than a single whole-symbol
   *  alias, if it has no current term yet, or if `ctx` doesn't apply
   *  to it (a COMBINATIONAL write to a non-wire symbol, or any
   *  INITIAL write, isn't needed by any currently-supported
   *  construct).
   *  @param sel  the dynamic element-select LHS expression
   *  @param rhs_expr the assignment's right-hand side
   *  @param ctx  which kind of block this assignment is in
   *  @param condition accumulated path condition (for if/case nesting)
   *  @param prefix the current hierarchical name prefix
   */
  /** @param base_expr what is being written into
   *  @param index_expr the runtime position, in elements for a
   *         element select and in bits for an indexed range select
   *  @param write_type the type of the slice being written, whose
   *         width says how much of the base the write covers
   *  @param scale_by_width multiply the index by that width, which
   *         an element select needs and a `+:`/`-:` does not, since
   *         its base already counts in bits
   *  @param pos_bias added to the position afterwards: zero for
   *         `+:`, and one less than the width, negated, for `-:`,
   *         whose base names the *top* of the range
   */
  void process_dynamic_write(const slang::ast::Expression & base_expr,
                             const slang::ast::Expression & index_expr,
                             const slang::ast::Type & write_type,
                             bool scale_by_width,
                             int64_t pos_bias,
                             const slang::ast::Expression & rhs_expr,
                             StmtContext ctx,
                             const smt::Term & condition,
                             const std::string & prefix);

  /** Build the array value a compile-time-constant assignment
   *  pattern denotes, or a null Term if `rhs_expr` is not one. A
   *  uniform fill stays a single constant array, with a Store only
   *  for each element that differs from the first. */
  smt::Term constant_array_term(
      const slang::ast::Expression & rhs_expr,
      const slang::ast::FixedSizeUnpackedArrayType & arr,
      const smt::Sort & array_sort);

  /** The same for an already-evaluated value, which is how it
   *  recurses through a multi-dimensional pattern. */
  smt::Term constant_array_value(
      const slang::ConstantValue & cv,
      const slang::ast::FixedSizeUnpackedArrayType & arr,
      const smt::Sort & array_sort);

  /** The value a whole-array target carries so far in the block being
   *  walked, or the array's own term if nothing has written it yet.
   *  A clocked block accumulates by term, the others by symbol. */
  smt::Term array_pending_value(const slang::ast::Symbol * sym,
                                const smt::Term & state_term,
                                StmtContext ctx);

  /** Record `value` as what that target now carries. */
  void record_array_pending(const slang::ast::Symbol * sym,
                            const smt::Term & state_term,
                            StmtContext ctx,
                            const smt::Term & value);

  /** Handle an assignment whose target is one element of an unpacked
   *  array, or a bit range inside one (`mem[i] <= v`, `mem[i][3:0] <=
   *  v`, `mem[i].f <= v`).  An element is not a bit range of its base,
   *  so none of these can go through commit_write(); each becomes a
   *  Store, over a Select-and-splice when the write is narrower than
   *  the element.
   *
   *  Returns false if `lhs_expr` does not name an unpacked-array
   *  element at all.  Any element write it cannot model throws:
   *  leaving an array unconstrained would read as an arbitrary value
   *  rather than as a missing feature. */
  bool process_array_element_assign(
      const slang::ast::Expression & lhs_expr,
      const slang::ast::Expression & rhs_expr,
      StmtContext ctx,
      const smt::Term & condition,
      const std::string & prefix,
      const smt::Term & rhs_override = smt::Term());

  /** Handle an assignment whose target is a whole unpacked array
   *  (`mem <= '0`), which is neither a bit range nor an element and so
   *  cannot go through resolve_lvalue()/LValueDesc.
   *
   *  Only a compile-time-constant right-hand side is supported: it
   *  becomes a constant array, with a Store for each element that
   *  differs from the first, so a uniform fill stays a single term.
   *  @return true if the assignment was handled here; false if it is
   *          not a whole-array assignment at all, leaving the caller's
   *          ordinary paths to deal with it
   */
  /** An assignment pattern whose elements are not all
   *  elaboration-time constants: one term per element, stored over
   *  `seed`. Every element is written, so `seed` only supplies the
   *  sort -- pass the array's own term. A constant array cannot
   *  serve here even for a uniform pattern, since not every solver
   *  accepts one built from a symbolic value.
   *
   *  Returns a null Term if `expr` is not an assignment pattern of
   *  the right shape. */
  smt::Term array_from_pattern(
      const slang::ast::Expression & expr,
      const slang::ast::FixedSizeUnpackedArrayType & arr,
      const smt::Term & seed,
      const std::string & prefix);

  bool process_whole_array_assign(const slang::ast::Expression & lhs_expr,
                                  const slang::ast::Expression & rhs_expr,
                                  StmtContext ctx,
                                  const smt::Term & condition,
                                  const std::string & prefix);

  /** Re-derive `symbol_table_.loop_var_terms()[&sym]` from `sym`'s
   *  current constant value in expr_encoder_.eval_ctx() (after a
   *  for-loop step, a while/repeat/foreach iteration, or a plain
   *  assignment to a compile-time-unrolled local).  Throws if `sym`
   *  isn't a currently-bound integer local. */
  void refresh_loop_var_term(const slang::ast::ValueSymbol & sym);

  /** The return variable of the subroutine being inlined, which a
   *  `return` binds. Null outside one, and also inside a task, which
   *  has no return value. */
  const slang::ast::Symbol * current_return_var_ = nullptr;

  /** Whether a subroutine body is being inlined at all, which is
   *  what decides that a write reaching outside it is an error --
   *  a task has no return variable, so that cannot be the test. */
  bool in_subroutine_ = false;

  /** Tasks currently being inlined, so recursion is reported rather
   *  than expanded forever. Functions are tracked separately, by
   *  ExprEncoder, which owns their call sites. */
  std::unordered_set<const slang::ast::Symbol *> inlining_tasks_;

  SymbolTable & symbol_table_;
  ExprEncoder & expr_encoder_;
  AssertionWalker & assertion_walker_;
  FunctionalTransitionSystem & fts_;
  const smt::SmtSolver & solver_;
};

}  // namespace pono
