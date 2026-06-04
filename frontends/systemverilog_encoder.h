/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   None
 ** This file is part of the pono project.
 ** Copyright (c) 2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief SystemVerilog frontend encoder using the slang library.
 **
 ** Parses and elaborates SystemVerilog designs via slang, then converts
 ** the elaborated representation into Pono's FunctionalTransitionSystem
 ** and properties.
 **
 **/

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "core/fts.h"
#include "smt-switch/smt.h"

// Forward declarations for slang types to avoid exposing slang headers
// to all translation units.
namespace slang::ast {
class Compilation;
class Type;
class Expression;
class Statement;
class Symbol;
class Scope;
class InstanceSymbol;
class InstanceBodySymbol;
class ProceduralBlockSymbol;
class ContinuousAssignSymbol;
class PortSymbol;
class EvalContext;
}  // namespace slang::ast

namespace pono {

class SystemVerilogEncoder
{
 public:
  /** Construct an SystemVerilogEncoder that parses a SystemVerilog file and populates
   *  the given FunctionalTransitionSystem.
   *
   *  Supported SystemVerilog constructs (synthesizable subset):
   *    - Module port declarations (input/output/inout)
   *    - always_ff blocks with non-blocking assignments (state elements)
   *    - always_comb blocks and continuous assignments (combinational logic)
   *    - initial blocks (initial state constraints)
   *    - Basic SVA immediate assertions (converted to properties)
   *    - Bitvector types (logic, bit, reg with packed dimensions)
   *
   *  @param filename path to the SystemVerilog source file
   *  @param fts the transition system to populate
   */
  SystemVerilogEncoder(std::string filename, FunctionalTransitionSystem & fts);

  ~SystemVerilogEncoder();

  /** @return the vector of properties (negated assertions) found in the design
   */
  const smt::TermVec & propvec() const { return propvec_; }

 private:
  // ---------- Encoding pipeline ----------

  /** Top-level encoding: parse file, elaborate, and walk the design. */
  void encode(const std::string & filename);

  /** Process a top-level module instance. */
  void process_module(const slang::ast::InstanceSymbol & inst);

  /** First pass: declare all variables (state vars, inputs, wires).
   *  Walks ports and internal variable declarations.
   */
  void declare_variables(const slang::ast::InstanceBodySymbol & body);

  /** Declare just the internal (non-port) variables of `body`.  Used
   *  when descending into a child instance, whose ports have already
   *  been bound through the port-connection map.
   */
  void declare_variables_internal(
      const slang::ast::InstanceBodySymbol & body);

  /** Declare a single port as an input or output variable. */
  void process_port(const slang::ast::PortSymbol & port);

  /** Second pass: process all behavioral and structural assignments.
   *  Walks always blocks, continuous assigns, and initial blocks.
   */
  void process_assignments(const slang::ast::InstanceBodySymbol & body);

  /** Process an always_ff block to extract next-state update functions. */
  void process_always_ff(const slang::ast::ProceduralBlockSymbol & proc);

  /** Process an always_comb block to extract combinational definitions. */
  void process_always_comb(const slang::ast::ProceduralBlockSymbol & proc);

  /** Process an initial block to extract initial-state constraints. */
  void process_initial(const slang::ast::ProceduralBlockSymbol & proc);

  /** Process a continuous assignment (assign statement). */
  void process_continuous_assign(
      const slang::ast::ContinuousAssignSymbol & ca);

  // ---------- Statement processing ----------

  /** Context for statement processing: whether we are building next-state
   *  updates (always_ff), combinational definitions (always_comb), or
   *  initial constraints.
   */
  enum class StmtContext
  {
    NEXT_STATE,     ///< Inside always_ff: build next-state functions
    COMBINATIONAL,  ///< Inside always_comb: build combinational definitions
    INITIAL         ///< Inside initial: build init constraints
  };

  /** Recursively process a statement, extracting assignments.
   *  @param stmt the slang statement to process
   *  @param ctx  what kind of block we are in
   *  @param condition accumulated path condition (for if/case nesting)
   */
  void process_statement(const slang::ast::Statement & stmt,
                         StmtContext ctx,
                         const smt::Term & condition);

  // ---------- Expression conversion ----------

  /** Convert a slang expression to an SMT term.
   *  Handles operators, literals, variable references, concatenation,
   *  bit-selects, ternary, etc.
   *  @param expr the slang expression
   *  @return the corresponding SMT term
   */
  smt::Term expr_to_term(const slang::ast::Expression & expr);

  // ---------- Type conversion ----------

  /** Convert a slang type to an SMT sort.
   *  @param type the slang type
   *  @return the corresponding SMT sort (BV or BOOL)
   */
  smt::Sort type_to_sort(const slang::ast::Type & type);

  // ---------- Helpers ----------

  /** Look up the SMT term for a slang symbol.
   *  @param sym pointer to the slang symbol
   *  @return the SMT term, or throws if not found
   */
  smt::Term lookup_symbol(const slang::ast::Symbol * sym) const;

  /** Build the hierarchical name for a symbol.
   *  @param name the local name
   *  @return the fully-qualified name with prefix
   */
  std::string make_name(const std::string & name) const;

  /** Ensure a term has the expected bit-width, zero-extending or
   *  truncating as needed.
   *  @param t the term to resize
   *  @param target_width the desired width
   *  @return the resized term
   */
  smt::Term resize_to(const smt::Term & t, uint64_t target_width);

  /** Build a partial-write term: take the full-width `base` and
   *  return a term equal to `base` everywhere except bits
   *  [lo .. hi], which take their values from `slice`.  A full-width
   *  write (lo == 0 && hi == width(base)-1) is just `slice`.
   */
  smt::Term replace_bits(const smt::Term & base,
                         const smt::Term & slice,
                         uint64_t lo, uint64_t hi);

  /** Pre-scan an always_ff body to identify non-blocking assignment
   *  targets as state variable symbols.
   *  @param body the statement body of the always_ff block
   */
  void pre_scan_always_ff(const slang::ast::Statement & body);

  /** Pre-scan a combinational always_comb body to identify blocking
   *  assignment targets as combinational wire symbols.
   *  @param body the statement body of the always_comb block
   */
  void pre_scan_always_comb(const slang::ast::Statement & body);

  /** Pre-scan a child instance to identify any parent-side variables
   *  that are driven by the child's output ports; those become wires
   *  in the parent's transition system.  Also recurses into nested
   *  instances.
   *  @param inst the child instance to scan
   */
  void pre_scan_instance(const slang::ast::InstanceSymbol & inst);

  /** Process a child module instance: register its port connections,
   *  then walk its body (continuous assigns, always blocks, nested
   *  instances) so the parent's transition system inherits all of
   *  the child's behavior.
   *  @param inst the child instance to process
   */
  void process_instance(const slang::ast::InstanceSymbol & inst);

  /** Invoke `fn` for every concrete member of `scope`, recursing
   *  through generate-block scopes so the caller sees the unrolled
   *  per-iteration members directly.  Uninstantiated generate
   *  blocks (the unselected arm of a generate-if / generate-case)
   *  are skipped.
   *  @param scope the scope whose members should be visited
   *  @param fn   callback invoked once per concrete member symbol
   */
  template <typename Fn>
  void walk_members(const slang::ast::Scope & scope, Fn && fn);

  /** Lazily construct and return the encoder's slang EvalContext.
   *  Used to evaluate compile-time-constant loop bounds and step
   *  expressions when unrolling procedural `for` loops.
   */
  slang::ast::EvalContext & eval_ctx();

  // ---------- Data members ----------

  FunctionalTransitionSystem & fts_;
  const smt::SmtSolver & solver_;

  // Unique pointer to slang compilation (hidden from header users)
  std::unique_ptr<slang::ast::Compilation> compilation_;

  // Map from slang symbol pointer to the corresponding SMT term.
  // For state variables, this maps to the current-state term.
  // For wires (continuous-assign / always_comb targets), this maps to
  // the defining SMT expression (macro substitution).
  std::unordered_map<const slang::ast::Symbol *, smt::Term> symbol_to_term_;

  // Set of symbols that are state variables (assigned in always_ff).
  // Populated during the first pass so that the second pass knows
  // which variables are registers vs. wires.
  std::unordered_set<const slang::ast::Symbol *> state_var_symbols_;

  // Set of symbols that are combinational wires (assigned in
  // continuous_assign or always_comb).  These are macro-substituted:
  // their term in symbol_to_term_ is the defining expression itself,
  // not a fresh state or input variable.
  std::unordered_set<const slang::ast::Symbol *> wire_symbols_;

  // Map from a child instance's port internal symbol to the parent-side
  // variable symbol that the port connects to.  Populated while
  // processing a child instance and consulted when resolving the LHS
  // of a continuous-assign or always_comb statement so writes to a
  // child output port redirect to the parent-side wire.
  std::unordered_map<const slang::ast::Symbol *, const slang::ast::Symbol *>
      port_output_aliases_;

  // Symbols in pending_comb_updates_ that came from an alias-redirect
  // through port_output_aliases_; their term must be named in
  // parent_prefix_ rather than prefix_ when the always_comb block
  // commits.
  std::unordered_set<const slang::ast::Symbol *> pending_comb_aliased_;

  // For always_ff processing: accumulated conditional next-state updates.
  // Maps state variable term -> conditional next-state expression.
  // After processing the block, these are committed via assign_next().
  std::unordered_map<smt::Term, smt::Term> pending_next_updates_;

  // For always_comb processing: accumulated wire-define expressions
  // keyed by the LHS wire symbol.  After processing the block, the
  // accumulated term is stored in symbol_to_term_.
  std::unordered_map<const slang::ast::Symbol *, smt::Term>
      pending_comb_updates_;

  // Properties extracted from SVA assert statements.
  smt::TermVec propvec_;

  // Hierarchical name prefix for the current module.
  std::string prefix_;

  // Hierarchical name prefix for the *parent* of the module currently
  // being processed (i.e., the prefix of the scope where output-port
  // aliases live).  Updated as we recurse into / out of child
  // instances so wires redirected via port_output_aliases_ get named
  // in their owning scope rather than the child's.
  std::string parent_prefix_;

  // Procedural-loop counter bindings: when unrolling a `for` loop,
  // each loop variable's slang VariableSymbol is mapped to a BV
  // term representing its current iteration value.  Consulted by
  // lookup_symbol so the loop body's references resolve to the
  // right constant in each unrolled iteration.
  std::unordered_map<const slang::ast::Symbol *, smt::Term>
      loop_var_terms_;

  // Slang evaluation context, lazily constructed on first use to
  // evaluate constant bounds and step expressions during for-loop
  // unrolling.  Owned via unique_ptr because EvalContext is not
  // default-constructible and is only forward-declared here.
  std::unique_ptr<slang::ast::EvalContext> eval_ctx_;

  // Stashed "current value of the LHS" used when expanding compound
  // assignments (`x &= y`, `x += y`, ...).  Slang represents the
  // implicit self-reference in the RHS as an
  // ExpressionKind::LValueReference; expr_to_term returns this term
  // for that case.  Set just before converting a compound RHS and
  // cleared right after, with save/restore for nested contexts.
  smt::Term current_lvalue_term_;
};

}  // namespace pono
