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
#include <vector>

#include "core/fts.h"
#include "smt-switch/smt.h"

namespace slang::ast {
class InstanceBodySymbol;
class PortSymbol;
class Scope;
class VariableSymbol;
}  // namespace slang::ast

namespace pono {

class SymbolTable;
struct ResolvedAliasPiece;

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
   *  been bound through the port-connection map -- and when descending
   *  into a checker instance, whose formal ports need no declaration
   *  at all (slang's own elaboration substitutes them). Takes a plain
   *  Scope (not specifically an InstanceBodySymbol) so it can also
   *  drive a CheckerInstanceBodySymbol's own local Variable/Net
   *  members, since walk_members() needs nothing InstanceBodySymbol-
   *  specific either.
   *  @param body the instance (or checker-instance) body to declare
   *         internal variables for
   *  @param prefix the hierarchical name prefix for `body`
   */
  void declare_variables_internal(const slang::ast::Scope & body,
                                  const std::string & prefix);

  /** Declare a single port as an input or output variable.
   *  @param port the port symbol to declare
   *  @param prefix the hierarchical name prefix for the enclosing instance
   */
  void process_port(const slang::ast::PortSymbol & port,
                    const std::string & prefix);

 private:
  /** Give an output-port-aliased register that covers only part of
   *  its target a state var of its own, and splice its bits into the
   *  target rather than aliasing the write onto it.
   *
   *  Aliasing cannot work once a target has several contributors:
   *  each sibling instance would claim the whole target's next-state
   *  function. Splicing composes instead, which is exactly how a comb
   *  wire driven from several sibling instances is already assembled
   *  -- including the ordering it inherits, where a read of the
   *  target before its contributors have been walked fails rather
   *  than seeing a half-built value.
   *
   *  Records what each contributor covers so the coverage check after
   *  every instance has been processed can tell whether the target
   *  ended up whole.
   */
  void splice_aliased_register(const slang::ast::VariableSymbol & var,
                               const std::vector<ResolvedAliasPiece> & pieces,
                               const std::string & prefix);

  SymbolTable & symbol_table_;
  FunctionalTransitionSystem & fts_;
  const smt::SmtSolver & solver_;
};

}  // namespace pono
