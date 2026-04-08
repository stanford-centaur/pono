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

#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "core/fts.h"
#include "smt-switch/smt.h"
#include "utils/exceptions.h"

// Forward declarations for slang types to avoid exposing slang headers
// to all translation units.
namespace slang::ast {
class Compilation;
class Type;
class Expression;
class Statement;
class Symbol;
class InstanceSymbol;
class InstanceBodySymbol;
class ProceduralBlockSymbol;
class ContinuousAssignSymbol;
class PortSymbol;
}  // namespace slang::ast

namespace pono {

class SVEncoder
{
 public:
  /** Construct an SVEncoder that parses a SystemVerilog file and populates
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
  SVEncoder(std::string filename, FunctionalTransitionSystem & fts);

  ~SVEncoder();

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

  /** Pre-scan an always_ff body to identify non-blocking assignment
   *  targets as state variable symbols.
   *  @param body the statement body of the always_ff block
   */
  void pre_scan_always_ff(const slang::ast::Statement & body);

  // ---------- Data members ----------

  FunctionalTransitionSystem & fts_;
  const smt::SmtSolver & solver_;

  // Unique pointer to slang compilation (hidden from header users)
  std::unique_ptr<slang::ast::Compilation> compilation_;

  // Map from slang symbol pointer to the corresponding SMT term.
  // For state variables, this maps to the current-state term.
  std::unordered_map<const slang::ast::Symbol *, smt::Term> symbol_to_term_;

  // Set of symbols that are state variables (assigned in always_ff).
  // Populated during the first pass so that the second pass knows
  // which variables are registers vs. wires.
  std::unordered_set<const slang::ast::Symbol *> state_var_symbols_;

  // For always_ff processing: accumulated conditional next-state updates.
  // Maps state variable term -> conditional next-state expression.
  // After processing the block, these are committed via assign_next().
  std::unordered_map<smt::Term, smt::Term> pending_next_updates_;

  // Properties extracted from SVA assert statements.
  smt::TermVec propvec_;

  // Hierarchical name prefix for the current module.
  std::string prefix_;
};

}  // namespace pono
