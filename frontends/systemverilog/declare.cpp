/*! \file declare.cpp
 *  \brief SystemVerilogEncoder's variable-declaration pass: creates the
 *         SMT term (state var / input var) for each port and internal
 *         variable/net, guided by the pre-scan classification.
 */
#include "frontends/systemverilog/encoder.h"
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

void SystemVerilogEncoder::declare_variables(
    const slang::ast::InstanceBodySymbol & body)
{
  using namespace slang::ast;

  // Process ports first.
  for (auto port_sym : body.getPortList()) {
    if (port_sym->kind == SymbolKind::Port) {
      process_port(port_sym->as<PortSymbol>());
    }
  }

  declare_variables_internal(body);
}

void SystemVerilogEncoder::declare_variables_internal(
    const slang::ast::InstanceBodySymbol & body)
{
  using namespace slang::ast;

  // Process internal variable declarations (non-port variables).
  walk_members(body, [&](const Symbol & member) {
    if (member.kind == SymbolKind::Variable) {
      auto & var = member.as<VariableSymbol>();
      // Skip if already declared via port processing.
      if (symbol_to_term_.count(&var)) return;
      // Wires get their term assigned during combinational-assignment
      // processing (macro substitution), not declared upfront.
      if (wire_symbols_.count(&var)) return;
      // Output ports of a child instance: the port-internal Variable
      // appears here as a member of the child's body, but its term is
      // really the parent-side wire reached through the alias map --
      // skip declaring a separate term for it.  A register is the
      // exception: unlike a comb wire (whose term is filled in later
      // via macro substitution when its driving assignment is
      // processed), no later pass ever assigns a term to a bare
      // pass-through symbol, so a register's state var must be
      // created here, keyed under the fully resolved alias root.
      if (port_output_aliases_.count(&var)) {
        if (state_var_symbols_.count(&var)) {
          uint64_t var_w = var.getType().getBitWidth();
          auto pieces = resolve_output_alias_pieces(&var, 0, var_w - 1);
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
            if (!symbol_to_term_.count(piece.sym)) {
              const Symbol * root = piece.sym;
              // Matches the pre-existing (single-target) naming
              // exactly when there's only one piece; disambiguated by
              // the target's own name for a concatenation-target
              // connection, where reusing `var.name` for every piece
              // would collide.
              string name =
                  make_name(pieces.size() > 1
                                ? string(var.name) + "_" + string(root->name)
                                : string(var.name));
              Sort sort = type_to_sort(root->as<ValueSymbol>().getType());
              Term sv = fts_.make_statevar(name, sort);
              symbol_to_term_[root] = sv;
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

      string name = make_name(string(var.name));
      Sort sort = type_to_sort(var.getType());

      if (state_var_symbols_.count(&var)) {
        // This is a register: create a state variable.
        Term sv = fts_.make_statevar(name, sort);
        symbol_to_term_[&var] = sv;
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
        symbol_to_term_[&var] = iv;
        fts_.name_term(name, iv);
        logger.log(2,
                   "SystemVerilogEncoder: undriven var {} : bv{}",
                   name,
                   sort->get_width());
      }
    } else if (member.kind == SymbolKind::Net) {
      auto & net = member.as<NetSymbol>();
      if (symbol_to_term_.count(&net)) return;
      if (wire_symbols_.count(&net)) return;
      if (port_output_aliases_.count(&net)) return;

      string name = make_name(string(net.name));
      Sort sort = type_to_sort(net.getType());

      Term iv = fts_.make_inputvar(name, sort);
      symbol_to_term_[&net] = iv;
      fts_.name_term(name, iv);
      logger.log(
          2, "SystemVerilogEncoder: net {} : bv{}", name, sort->get_width());
    }
  });
}

void SystemVerilogEncoder::process_port(const slang::ast::PortSymbol & port)
{
  using namespace slang::ast;

  string name = make_name(string(port.name));
  Sort sort = type_to_sort(port.getType());

  const Symbol * internal = port.internalSymbol;
  if (!internal) {
    // Port with no internal symbol -- create based on direction.
    if (port.direction == ArgumentDirection::In) {
      Term iv = fts_.make_inputvar(name, sort);
      symbol_to_term_[&port] = iv;
      fts_.name_term(name, iv);
      logger.log(2,
                 "SystemVerilogEncoder: input port {} : bv{}",
                 name,
                 sort->get_width());
    } else {
      // Output/inout without internal symbol: treat as state or wire
      // depending on whether it was found in always_ff pre-scan.
      if (state_var_symbols_.count(&port)) {
        Term sv = fts_.make_statevar(name, sort);
        symbol_to_term_[&port] = sv;
        fts_.name_term(name, sv);
      } else {
        Term iv = fts_.make_inputvar(name, sort);
        symbol_to_term_[&port] = iv;
        fts_.name_term(name, iv);
      }
    }
    return;
  }

  // Port has an internal symbol -- use it.
  if (port.direction == ArgumentDirection::In) {
    Term iv = fts_.make_inputvar(name, sort);
    symbol_to_term_[internal] = iv;
    symbol_to_term_[&port] = iv;
    fts_.name_term(name, iv);
    logger.log(2,
               "SystemVerilogEncoder: input port {} : bv{}",
               name,
               sort->get_width());
  } else {
    // Output or inout: classify based on driver kind.
    if (state_var_symbols_.count(internal)) {
      Term sv = fts_.make_statevar(name, sort);
      symbol_to_term_[internal] = sv;
      symbol_to_term_[&port] = sv;
      fts_.name_term(name, sv);
      logger.log(2,
                 "SystemVerilogEncoder: output port (reg) {} : bv{}",
                 name,
                 sort->get_width());
    } else if (wire_symbols_.count(internal)) {
      // Combinational output port: defer term creation to
      // process_continuous_assign / process_always_comb, which will
      // populate symbol_to_term_ for `internal`.  Lookups via the port
      // symbol fall back to the internal symbol's entry.
      logger.log(
          2, "SystemVerilogEncoder: output port (wire) {}: deferred", name);
    } else {
      Term iv = fts_.make_inputvar(name, sort);
      symbol_to_term_[internal] = iv;
      symbol_to_term_[&port] = iv;
      fts_.name_term(name, iv);
      logger.log(2,
                 "SystemVerilogEncoder: output port (undriven) {} : bv{}",
                 name,
                 sort->get_width());
    }
  }
}

}  // namespace pono
