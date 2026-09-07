/*!
 * \file expr_encoder.h
 * \brief expr_to_term(): converts slang AST expressions to SMT terms.
 * \author Áron Ricardo Perez-Lopez
 * \date 2026
 * \copyright See the LICENSE file in the top-level source directory.
 *
 * ExprEncoder has no dependency on SystemVerilogEncoder itself -- only on
 * SymbolTable (to resolve a NamedValue/HierarchicalValue read) and Tableau
 * (for the sampled-value system functions, $past/$stable/$changed/$rose/
 * $fell). It's called extensively by whichever class walks procedural
 * statements (confirmed one-directional: nothing in this file ever calls
 * back into statement processing).
 */
#pragma once

#include <memory>
#include <string>

#include "smt-switch/smt.h"

namespace slang::ast {
class Compilation;
class EvalContext;
class Expression;
}  // namespace slang::ast

namespace pono {

class SymbolTable;
class Tableau;

class ExprEncoder
{
 public:
  ExprEncoder(SymbolTable & symbol_table,
              Tableau & tableau,
              const smt::SmtSolver & solver);
  ~ExprEncoder();

  /** Must be called once, after the slang Compilation exists (i.e. after
   *  parsing/elaboration), and before the first call to eval_ctx() or to
   *  an expr_to_term() case that constant-folds via eval_ctx().
   */
  void bind_compilation(slang::ast::Compilation & compilation);

  /** Convert a slang expression to an SMT term.
   *  Handles operators, literals, variable references, concatenation,
   *  bit-selects, ternary, etc.
   *  @param expr the slang expression
   *  @param prefix the caller's current hierarchical name prefix, used
   *         only for naming hidden latches introduced by the sampled-
   *         value system functions ($past/$stable/$changed/$rose/$fell)
   *  @return the corresponding SMT term
   */
  smt::Term expr_to_term(const slang::ast::Expression & expr,
                         const std::string & prefix);

  /** Sets the term expr_to_term()'s LValueReference case returns (the
   *  implicit self-reference in a compound assignment's RHS, e.g. `x +=
   *  y`), and returns the previous value so the caller can restore it
   *  once the compound RHS has been converted.
   */
  smt::Term set_current_lvalue_term(const smt::Term & t);

  /** Lazily construct and return the shared slang EvalContext, used
   *  throughout statement/expression processing to evaluate compile-
   *  time-constant expressions -- not just loop bounds and step
   *  expressions when unrolling procedural loops, but also things like
   *  case-statement selectors/patterns, constant array/bit-select
   *  indices, and other constant-foldable operands (e.g. a `**` exponent
   *  or `$past`'s cycle-count argument).
   */
  slang::ast::EvalContext & eval_ctx();

 private:
  SymbolTable & symbol_table_;
  Tableau & tableau_;
  const smt::SmtSolver & solver_;
  slang::ast::Compilation * compilation_ = nullptr;

  // Stashed "current value of the LHS" used when expanding compound
  // assignments (`x &= y`, `x += y`, ...).  Slang represents the
  // implicit self-reference in the RHS as an
  // ExpressionKind::LValueReference; expr_to_term returns this term for
  // that case.  Set via set_current_lvalue_term() just before
  // converting a compound RHS and restored right after, with
  // save/restore for nested contexts.
  smt::Term current_lvalue_term_;

  // Slang evaluation context, lazily constructed on first use to
  // evaluate constant bounds and step expressions during for-loop
  // unrolling.  Owned via unique_ptr because EvalContext is not
  // default-constructible and is only forward-declared here.
  std::unique_ptr<slang::ast::EvalContext> eval_ctx_;
};

}  // namespace pono
