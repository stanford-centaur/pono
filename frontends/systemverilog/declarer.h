/*!
 * \file declarer.h
 * \brief Creates SMT terms for ports, registers, and free variables.
 * \author Áron Ricardo Perez-Lopez
 * \date 2026
 * \copyright See the LICENSE file in the top-level source directory.
 *
 * Declarer owns the variable-declaration pass: ports are declared first
 * (process_port()), then remaining internal Variable and Net symbols
 * (declare_variables_internal()). Each symbol is classified using
 * SymbolTable's pre-scan results: registers become state vars, undriven
 * signals become free input vars, and combinational wires are skipped
 * here -- their terms are filled in later, during macro-substitution in
 * continuous-assignment/always_comb processing.
 *
 * Depends only on SymbolTable (for the pre-scan classifications and
 * symbol_to_term_) and the FunctionalTransitionSystem/solver it declares
 * terms into -- both independent of the rest of the encoder, so this
 * class is too: it holds no reference back to SystemVerilogEncoder. The
 * hierarchical name prefix is a parameter to every entry point rather
 * than ambient state, matching SymbolTable's own design.
 */
#pragma once

#include <string>

#include "core/fts.h"
#include "smt-switch/smt.h"

namespace slang::ast {
class InstanceBodySymbol;
class PortSymbol;
}  // namespace slang::ast

namespace pono {

class SymbolTable;

class Declarer
{
 public:
  Declarer(SymbolTable & symbol_table,
           FunctionalTransitionSystem & fts,
           const smt::SmtSolver & solver);

  /** First pass: declare state vars and inputs.  Wires are skipped --
   *  they get their term assigned later during combinational-assignment
   *  processing.  Walks ports and internal variable declarations.
   *  @param body the instance body to declare variables for
   *  @param prefix the hierarchical name prefix for `body`
   */
  void declare_variables(const slang::ast::InstanceBodySymbol & body,
                         const std::string & prefix);

  /** Declare just the internal (non-port) variables of `body`.  Used
   *  when descending into a child instance, whose ports have already
   *  been bound through the port-connection map.
   *  @param body the instance body to declare internal variables for
   *  @param prefix the hierarchical name prefix for `body`
   */
  void declare_variables_internal(const slang::ast::InstanceBodySymbol & body,
                                  const std::string & prefix);

  /** Declare a single port as an input or output variable.
   *  @param port the port symbol to declare
   *  @param prefix the hierarchical name prefix for the enclosing instance
   */
  void process_port(const slang::ast::PortSymbol & port,
                    const std::string & prefix);

 private:
  SymbolTable & symbol_table_;
  FunctionalTransitionSystem & fts_;
  const smt::SmtSolver & solver_;
};

}  // namespace pono
