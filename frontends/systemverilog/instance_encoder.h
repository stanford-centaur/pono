/*!
 * \file instance_encoder.h
 * \brief Per-instance pass: continuous assigns, procedural blocks, instances.
 * \author Áron Ricardo Perez-Lopez
 * \date 2026
 * \copyright See the LICENSE file in the top-level source directory.
 *
 * InstanceEncoder processes one module instance body in two passes so
 * combinational definitions are visible to their consumers: continuous
 * assigns and always_comb (including plain `always` blocks with no
 * nonblocking targets) first, then always_ff, always_latch, and initial
 * blocks second. Legacy `initial forever @(...)` is redirected to the
 * same NEXT_STATE handling as always_ff/always_latch. process_instance()
 * recurses into child instances, passing down a freshly computed
 * hierarchical name prefix (and the parent's own prefix, for output-port
 * aliases) rather than mutating shared state -- since these are plain
 * parameters rather than an ambient member, no save/restore is needed
 * around the recursion the way symbol_table_'s driver-prefix bookkeeping
 * still needs it. `checker`, `specify`, `program`, and `defparam`
 * constructs are recognized and skipped as simulation-only or
 * functionally-inert.
 *
 * Implements SymbolTable::DriverResolver for real (see symbol_table.h):
 * SymbolTable's on-demand wire resolution is a genuine, load-bearing
 * mutual dependency with this class's own driver processing, so this is
 * the one seam in the whole rearchitecture with an actual back-and-forth
 * dependency, funneled through that one small interface rather than a
 * back-reference to this whole class.
 *
 * Depends on SymbolTable, Declarer, StatementEncoder, and ExprEncoder --
 * all already independent of the rest of the encoder, so this class is
 * too: it holds no reference back to SystemVerilogEncoder. Needs its own
 * Compilation dependency (bound via bind_compilation(), the same two-
 * phase pattern ExprEncoder uses, since a Compilation does not exist yet
 * when this class is constructed) purely to resolve each instance body's
 * `default disable iff` via Compilation::getDefaultDisable() before
 * dispatching into StatementEncoder.
 */
#pragma once

#include <string>

#include "core/fts.h"
#include "frontends/systemverilog/symbol_table.h"
#include "smt-switch/smt.h"

namespace slang::ast {
class Compilation;
class ContinuousAssignSymbol;
class Expression;
class InstanceBodySymbol;
class InstanceSymbol;
class ProceduralBlockSymbol;
class Scope;
class Statement;
}  // namespace slang::ast

namespace pono {

class Declarer;
class ExprEncoder;
class StatementEncoder;

class InstanceEncoder : private SymbolTable::DriverResolver
{
 public:
  InstanceEncoder(SymbolTable & symbol_table,
                  Declarer & declarer,
                  StatementEncoder & statement_encoder,
                  ExprEncoder & expr_encoder,
                  FunctionalTransitionSystem & fts,
                  const smt::SmtSolver & solver);

  /** Bind the slang Compilation this instance encoder should resolve
   *  `default disable iff` expressions against. Must be called once,
   *  after the Compilation is parsed and before process_assignments()
   *  is first called (see SystemVerilogEncoder::run()).
   */
  void bind_compilation(slang::ast::Compilation & compilation);

  /** Second pass: process all behavioral and structural assignments.
   *  Walks always blocks, continuous assigns, and initial blocks.
   *  @param body the instance body to process
   *  @param prefix the hierarchical name prefix for `body`
   *  @param parent_prefix the hierarchical name prefix of `body`'s
   *         enclosing scope (where output-port aliases live)
   */
  void process_assignments(const slang::ast::InstanceBodySymbol & body,
                           const std::string & prefix,
                           const std::string & parent_prefix);

  /** Process a child module instance: register its port connections,
   *  then walk its body (continuous assigns, always blocks, nested
   *  instances) so the parent's transition system inherits all of the
   *  child's behavior.
   *  @param inst the child instance to process
   *  @param prefix the hierarchical name prefix of `inst`'s enclosing
   *         scope (the child's own prefix is computed from this plus
   *         `inst`'s name)
   *  @param parent_prefix the hierarchical name prefix that `prefix`
   *         itself was computed relative to (becomes the child's own
   *         parent_prefix)
   */
  void process_instance(const slang::ast::InstanceSymbol & inst,
                        const std::string & prefix,
                        const std::string & parent_prefix);

 private:
  void process_always_ff(const slang::ast::ProceduralBlockSymbol & proc,
                         const std::string & prefix);

  /** Shared implementation of process_always_ff(): process `body` with
   *  StatementEncoder::StmtContext::NEXT_STATE and commit every
   *  accumulated pending_next_updates_ entry via assign_next(). Blocking-
   *  vs. nonblocking-assignment syntax makes no difference to this
   *  encoding -- only the StmtContext does -- so this is also reused
   *  directly for `always_latch` (whose writes use blocking `=`, but
   *  need exactly the same "holds previous value when not written"
   *  register semantics) and for a `forever @(...) ...` loop recognized
   *  as a legacy structural spelling of `always_ff`.
   */
  void process_next_state_body(const slang::ast::Statement & body,
                               const std::string & prefix);

  void process_always_comb(const slang::ast::ProceduralBlockSymbol & proc,
                           const std::string & prefix,
                           const std::string & parent_prefix);

  void process_initial(const slang::ast::ProceduralBlockSymbol & proc,
                       const std::string & prefix);

  /** Process a continuous assignment (assign statement). */
  void process_continuous_assign(const slang::ast::ContinuousAssignSymbol & ca,
                                 const std::string & prefix,
                                 const std::string & parent_prefix);

  /** Process one write from a continuous assignment: `lhs_expr` is
   *  either the assignment's whole LHS (the common case), or one
   *  operand of a concatenation-target LHS (`assign {hi, lo} = ...`),
   *  in which case `rhs` is already sliced down to that operand's own
   *  width.
   *  @param lhs_expr the (non-concatenation) lvalue expression being
   *                  written
   *  @param rhs the value to write, resized to lhs_expr's width
   *  @param rhs_signed whether to sign-extend (rather than zero-
   *         extend) `rhs` if it needs to grow to fit lhs_expr's width
   *         -- the caller passes the *source* expression's own
   *         signedness (false for a concatenation-target slice, which
   *         is positional bit-splicing, not a numeric value)
   */
  void process_continuous_assign_operand(
      const slang::ast::Expression & lhs_expr,
      const smt::Term & rhs,
      bool rhs_signed,
      const std::string & prefix,
      const std::string & parent_prefix);

  /** Process an always_comb block, or a continuous assign, unless it
   *  has already been processed (either by the normal walk reaching
   *  it, or because some earlier read forced it via
   *  resolve_wire_on_demand()).  Every call site that would otherwise
   *  call process_always_comb() / process_continuous_assign() directly
   *  goes through these instead, so a driver is never processed twice.
   */
  void process_always_comb_once(const slang::ast::ProceduralBlockSymbol & proc,
                                const std::string & prefix,
                                const std::string & parent_prefix);
  void process_continuous_assign_once(
      const slang::ast::ContinuousAssignSymbol & ca,
      const std::string & prefix,
      const std::string & parent_prefix);

  // ---------- SymbolTable::DriverResolver overrides ----------
  void resolve_continuous_assign(const slang::ast::ContinuousAssignSymbol & ca,
                                 const std::string & prefix,
                                 const std::string & parent_prefix) override;
  void resolve_always_comb(const slang::ast::ProceduralBlockSymbol & proc,
                           const std::string & prefix,
                           const std::string & parent_prefix) override;

  SymbolTable & symbol_table_;
  Declarer & declarer_;
  StatementEncoder & statement_encoder_;
  ExprEncoder & expr_encoder_;
  FunctionalTransitionSystem & fts_;
  const smt::SmtSolver & solver_;

  // Bound once by bind_compilation(); used only to resolve each
  // instance body's `default disable iff` via getDefaultDisable().
  slang::ast::Compilation * compilation_ = nullptr;

  // Scope of the module instance body currently being processed by
  // process_assignments().  Used to look up that module's `default
  // disable iff` condition (Compilation::getDefaultDisable walks up
  // the scope's parent chain).  Saved/restored around recursion into
  // child instances; note that process_instance() does not itself
  // update current_scope_ when it recurses into a child instance, so
  // assertions inside a child instance's own procedural blocks are
  // still resolved against the outer scope rather than the child's.
  const slang::ast::Scope * current_scope_ = nullptr;
};

}  // namespace pono
