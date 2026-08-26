/*!
 * \file declarer.cpp
 * \brief Creates SMT terms for ports, registers, and free variables.
 * \author Áron Ricardo Perez-Lopez
 * \date 2026
 * \copyright See the LICENSE file in the top-level source directory.
 *
 * A special case handles registers reached through an output-port alias
 * chain, splicing per-piece state vars keyed by the fully resolved alias
 * root; only full-width aliasing is supported today, partial-width
 * aliasing throws.
 */
#include "frontends/systemverilog/declarer.h"

#include <string>

#include "frontends/systemverilog/ast_helpers.h"
#include "frontends/systemverilog/bit_utils.h"
#include "frontends/systemverilog/symbol_table.h"
#include "slang/ast/Symbol.h"
#include "slang/ast/symbols/InstanceSymbols.h"
#include "slang/ast/symbols/PortSymbols.h"
#include "slang/ast/symbols/VariableSymbols.h"
#include "slang/ast/types/Type.h"
#include "smt-switch/smt.h"
#include "utils/exceptions.h"
#include "utils/logger.h"

using namespace smt;
using namespace std;

namespace pono {

Declarer::Declarer(SymbolTable & symbol_table,
                   FunctionalTransitionSystem & fts,
                   const smt::SmtSolver & solver)
    : symbol_table_(symbol_table), fts_(fts), solver_(solver)
{
}

void Declarer::declare_variables(const slang::ast::InstanceBodySymbol & body,
                                 const string & prefix)
{
  using namespace slang::ast;

  // Process ports first.
  for (auto port_sym : body.getPortList()) {
    if (port_sym->kind == SymbolKind::Port) {
      process_port(port_sym->as<PortSymbol>(), prefix);
    }
  }

  declare_variables_internal(body, prefix);
}

void Declarer::declare_variables_internal(
    const slang::ast::InstanceBodySymbol & body, const string & prefix)
{
  using namespace slang::ast;

  // Process internal variable declarations (non-port variables).
  // walk_members() takes the prefix by mutable reference (updating it
  // while descending into generate-for/instance-array child scopes,
  // then restoring it), so it needs its own local copy rather than
  // `prefix` itself.
  string walk_prefix = prefix;
  walk_members(body, walk_prefix, [&](const Symbol & member) {
    if (member.kind == SymbolKind::Variable) {
      auto & var = member.as<VariableSymbol>();
      // Skip if already declared via port processing.
      if (symbol_table_.symbol_to_term().count(&var)) return;
      // Wires get their term assigned during combinational-assignment
      // processing (macro substitution), not declared upfront.
      if (symbol_table_.wire_symbols().count(&var)) return;
      // Output ports of a child instance: the port-internal Variable
      // appears here as a member of the child's body, but its term is
      // really the parent-side wire reached through the alias map --
      // skip declaring a separate term for it.  A register is the
      // exception: unlike a comb wire (whose term is filled in later
      // via macro substitution when its driving assignment is
      // processed), no later pass ever assigns a term to a bare
      // pass-through symbol, so a register's state var must be
      // created here, keyed under the fully resolved alias root.
      if (symbol_table_.port_output_aliases().count(&var)) {
        if (symbol_table_.state_var_symbols().count(&var)) {
          uint64_t var_w = var.getType().getBitWidth();
          auto pieces =
              symbol_table_.resolve_output_alias_pieces(&var, 0, var_w - 1);
          for (auto & piece : pieces) {
            uint64_t target_w =
                piece.sym->as<ValueSymbol>().getType().getBitWidth();
            bool piece_full =
                (piece.target_lo == 0 && piece.target_hi + 1 == target_w);
            // A piece that only covers *part* of its own target's
            // width (e.g. one element of an instance array wired to a
            // slice of a shared register) has no splicing logic here
            // the way process_continuous_assign() has for a wire --
            // throw rather than silently never create a state var for
            // the shared target at all.
            if (!piece_full) {
              throw PonoException(
                  "SystemVerilogEncoder: register '" + string(var.name)
                  + "' is output-port-aliased to only part of '"
                  + string(piece.sym->name)
                  + "' (e.g. one element of an instance array wired to "
                    "a slice of a shared register) -- not supported");
            }
            if (!symbol_table_.symbol_to_term().count(piece.sym)) {
              const Symbol * root = piece.sym;
              // Matches the pre-existing (single-target) naming
              // exactly when there's only one piece; disambiguated by
              // the target's own name for a concatenation-target
              // connection, where reusing `var.name` for every piece
              // would collide.
              string name = symbol_table_.make_name(
                  walk_prefix,
                  pieces.size() > 1
                      ? string(var.name) + "_" + string(root->name)
                      : string(var.name));
              Sort sort =
                  type_to_sort(solver_, root->as<ValueSymbol>().getType());
              Term sv = fts_.make_statevar(name, sort);
              symbol_table_.symbol_to_term()[root] = sv;
              fts_.name_term(name, sv);
              logger.log(2,
                         "SystemVerilogEncoder: state var (aliased) {} : bv{}",
                         name,
                         sort->get_width());
            }
          }
        }
        return;
      }

      string name = symbol_table_.make_name(walk_prefix, string(var.name));
      Sort sort = type_to_sort(solver_, var.getType());

      if (symbol_table_.state_var_symbols().count(&var)) {
        // This is a register: create a state variable.
        Term sv = fts_.make_statevar(name, sort);
        symbol_table_.symbol_to_term()[&var] = sv;
        fts_.name_term(name, sv);
        logger.log(2,
                   "SystemVerilogEncoder: state var {} : bv{}",
                   name,
                   sort->get_width());
      } else {
        // No assignment found for this variable -- treat as a free
        // input.  This matches Verilog's "open" semantics where an
        // undriven net is unconstrained.
        Term iv = fts_.make_inputvar(name, sort);
        symbol_table_.symbol_to_term()[&var] = iv;
        fts_.name_term(name, iv);
        logger.log(2,
                   "SystemVerilogEncoder: undriven var {} : bv{}",
                   name,
                   sort->get_width());
      }
    } else if (member.kind == SymbolKind::Net) {
      auto & net = member.as<NetSymbol>();
      if (symbol_table_.symbol_to_term().count(&net)) return;
      if (symbol_table_.wire_symbols().count(&net)) return;
      if (symbol_table_.port_output_aliases().count(&net)) return;

      string name = symbol_table_.make_name(walk_prefix, string(net.name));
      Sort sort = type_to_sort(solver_, net.getType());

      Term iv = fts_.make_inputvar(name, sort);
      symbol_table_.symbol_to_term()[&net] = iv;
      fts_.name_term(name, iv);
      logger.log(
          2, "SystemVerilogEncoder: net {} : bv{}", name, sort->get_width());
    }
  });
}

void Declarer::process_port(const slang::ast::PortSymbol & port,
                            const string & prefix)
{
  using namespace slang::ast;

  string name = symbol_table_.make_name(prefix, string(port.name));
  Sort sort = type_to_sort(solver_, port.getType());

  const Symbol * internal = port.internalSymbol;
  if (!internal) {
    // Port with no internal symbol -- create based on direction.
    if (port.direction == ArgumentDirection::In) {
      Term iv = fts_.make_inputvar(name, sort);
      symbol_table_.symbol_to_term()[&port] = iv;
      fts_.name_term(name, iv);
      logger.log(2,
                 "SystemVerilogEncoder: input port {} : bv{}",
                 name,
                 sort->get_width());
    } else {
      // Output/inout without internal symbol: treat as a state var if
      // pre-scan found it driven by always_ff/always/always_latch,
      // otherwise as an (unconstrained) input var.
      if (symbol_table_.state_var_symbols().count(&port)) {
        Term sv = fts_.make_statevar(name, sort);
        symbol_table_.symbol_to_term()[&port] = sv;
        fts_.name_term(name, sv);
      } else {
        Term iv = fts_.make_inputvar(name, sort);
        symbol_table_.symbol_to_term()[&port] = iv;
        fts_.name_term(name, iv);
      }
    }
    return;
  }

  // Port has an internal symbol -- use it.
  if (port.direction == ArgumentDirection::In) {
    Term iv = fts_.make_inputvar(name, sort);
    symbol_table_.symbol_to_term()[internal] = iv;
    symbol_table_.symbol_to_term()[&port] = iv;
    fts_.name_term(name, iv);
    logger.log(2,
               "SystemVerilogEncoder: input port {} : bv{}",
               name,
               sort->get_width());
  } else {
    // Output or inout: classify based on driver kind.
    if (symbol_table_.state_var_symbols().count(internal)) {
      Term sv = fts_.make_statevar(name, sort);
      symbol_table_.symbol_to_term()[internal] = sv;
      symbol_table_.symbol_to_term()[&port] = sv;
      fts_.name_term(name, sv);
      logger.log(2,
                 "SystemVerilogEncoder: output port (reg) {} : bv{}",
                 name,
                 sort->get_width());
    } else if (symbol_table_.wire_symbols().count(internal)) {
      // Combinational output port: defer term creation to
      // process_continuous_assign / process_always_comb, which will
      // populate symbol_to_term_ for `internal` once its driving
      // assignment is processed.  Expressions elsewhere reference
      // `internal` directly (per slang, the port symbol itself is
      // never referenceable), so no entry for `&port` is needed.
      logger.log(
          2, "SystemVerilogEncoder: output port (wire) {}: deferred", name);
    } else {
      Term iv = fts_.make_inputvar(name, sort);
      symbol_table_.symbol_to_term()[internal] = iv;
      symbol_table_.symbol_to_term()[&port] = iv;
      fts_.name_term(name, iv);
      logger.log(2,
                 "SystemVerilogEncoder: output port (undriven) {} : bv{}",
                 name,
                 sort->get_width());
    }
  }
}

}  // namespace pono
