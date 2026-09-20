/*!
 * \file expr_encoder.h
 * \brief expr_to_term()/expr_to_bool(): slang AST expressions to SMT terms.
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

  /** Convert a slang expression to a bit-vector SMT term.
   *  Handles operators, literals, variable references, concatenation,
   *  bit-selects, ternary, etc.  Always BV-sorted, even for a 1-bit SV
   *  expression, since any expression can appear where a bit-vector is
   *  required (inside a concatenation, a part-select, a width-changing
   *  conversion).  Callers that only want the expression's truth value
   *  should prefer expr_to_bool().
   *  @param expr the slang expression
   *  @param prefix the caller's current hierarchical name prefix, used
   *         only for naming hidden latches introduced by the sampled-
   *         value system functions ($past/$stable/$changed/$rose/$fell)
   *  @return the corresponding SMT term, BV-sorted
   */
  smt::Term expr_to_term(const slang::ast::Expression & expr,
                         const std::string & prefix);

  /** Convert a slang expression to a Bool SMT term giving its SV truth
   *  value.  Equivalent to `expr_to_term(expr) != 0`, and exactly that
   *  for an expression whose value is a bit-vector -- but an expression
   *  that is already a predicate (a comparison, `&&`/`||`/`!`, an
   *  and/or reduction, $rose/$fell/$stable/$changed/$onehot/...) is
   *  returned directly, instead of being materialised as a 1-bit
   *  bit-vector only for the caller to compare it against zero again.
   *  @param expr the slang expression
   *  @param prefix see expr_to_term()
   *  @return the corresponding SMT term, Bool-sorted
   */
  smt::Term expr_to_bool(const slang::ast::Expression & expr,
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
  /** The single expression-conversion switch, returning `expr`'s value
   *  in its *natural* SMT sort: Bool for the expressions that really
   *  are predicates, a bit-vector for everything else.  Only the two
   *  public wrappers above call this; they adapt whichever sort comes
   *  back, which is what keeps the operator dispatch in one place and
   *  makes expr_to_bool() and expr_to_term() agree by construction.
   */
  /** Extract bits [lo, hi] of `val`, which may name bits the vector
   *  does not have. The LRM reads those as X, so they come back as a
   *  fresh unconstrained value rather than reaching the solver as an
   *  out-of-bounds Extract (which aborts) or being shifted in as
   *  zeros (which is silently wrong). */
  smt::Term extract_maybe_out_of_range(const smt::Term & val,
                                       uint64_t hi,
                                       uint64_t lo);

  smt::Term expr_to_term_or_bool(const slang::ast::Expression & expr,
                                 const std::string & prefix);

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
