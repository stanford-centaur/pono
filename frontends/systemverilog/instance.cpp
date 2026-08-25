/*!
 * \file instance.cpp
 * \brief Per-instance pass: continuous assigns, procedural blocks, instances.
 * \author Áron Ricardo Perez-Lopez
 * \date 2026
 * \copyright See the LICENSE file in the top-level source directory.
 *
 * Processes one module instance body in two passes so combinational definitions
 * are visible to their consumers: continuous assigns and always_comb (including
 * plain `always` blocks with no nonblocking targets) first, then always_ff,
 * always_latch, and initial blocks second. Legacy `initial forever @(...)` is
 * redirected to the same NEXT_STATE handling as always_ff/always_latch.
 * process_instance() recurses into child instances, rebinding hierarchical name
 * prefixes and port connections (including concatenation- and
 * array-element-target ports) as parent-side term reads (inputs) or
 * write-redirecting aliases (outputs), then undoes those bindings on return
 * since slang may share one InstanceBody across several instantiations.
 * `checker`, `specify`, `program`, and `defparam` constructs are recognized and
 * skipped as simulation-only or functionally-inert.
 */
#include <unordered_set>
#include <vector>

#include "frontends/systemverilog/ast_helpers.h"
#include "frontends/systemverilog/bit_utils.h"
#include "frontends/systemverilog/encoder.h"
#include "slang/ast/Compilation.h"
#include "slang/ast/EvalContext.h"
#include "slang/ast/Expression.h"
#include "slang/ast/SemanticFacts.h"
#include "slang/ast/Statement.h"
#include "slang/ast/Symbol.h"
#include "slang/ast/expressions/AssignmentExpressions.h"
#include "slang/ast/expressions/OperatorExpressions.h"
#include "slang/ast/symbols/BlockSymbols.h"
#include "slang/ast/symbols/CompilationUnitSymbols.h"
#include "slang/ast/symbols/InstanceSymbols.h"
#include "slang/ast/symbols/MemberSymbols.h"
#include "slang/ast/symbols/PortSymbols.h"
#include "slang/ast/types/Type.h"
#include "smt-switch/smt.h"
#include "utils/exceptions.h"
#include "utils/logger.h"

using namespace smt;
using namespace std;

namespace pono {

void SystemVerilogEncoder::process_assignments(
    const slang::ast::InstanceBodySymbol & body)
{
  using namespace slang::ast;

  // Track which module's scope we're in so assertions processed below
  // (in this function's own two walks) can look up that module's
  // `default disable iff` (if any).  Saved and restored around this
  // call; note that process_instance() below does not itself update
  // current_scope_ when it recurses into a child instance, so
  // assertions inside a child instance's own procedural blocks are
  // still resolved against this (outer) scope rather than the
  // child's.
  const Scope * saved_scope = current_scope_;
  current_scope_ = &body;

  // Combinational definitions are processed first so that wires have a
  // term assigned in symbol_to_term_ before any always_ff or initial
  // block (or assertion) tries to reference them.  Child instances are
  // walked here too -- a child's continuous assigns / always_comb
  // blocks may drive parent-side wires that downstream parent code
  // references.
  walk_members(body, prefix_, [&](const Symbol & member) {
    if (member.kind == SymbolKind::ContinuousAssign) {
      process_continuous_assign_once(member.as<ContinuousAssignSymbol>());
    } else if (member.kind == SymbolKind::ProceduralBlock) {
      auto & proc = member.as<ProceduralBlockSymbol>();
      if (proc.procedureKind == ProceduralBlockKind::AlwaysComb) {
        process_always_comb_once(proc);
      } else if (proc.procedureKind == ProceduralBlockKind::Always) {
        std::unordered_set<const Symbol *> targets;
        collect_nonblocking_targets(proc.getBody(), targets);
        if (targets.empty()) {
          process_always_comb_once(proc);
        }
      }
    } else if (member.kind == SymbolKind::Instance) {
      process_instance(member.as<InstanceSymbol>());
    } else if (member.kind == SymbolKind::CheckerInstance) {
      // `checker ... endchecker`, instantiated like a module: a
      // verification-only construct (SVA sequence/property checking
      // outside the DUT's own functional logic), so this instance and
      // everything inside it -- including any of its own assertions --
      // is simulation-only and has no functional-logic counterpart to
      // encode.
      logger.log(1,
                 "SystemVerilogEncoder: ignoring checker instance '{}' "
                 "(simulation-only construct)",
                 string(member.name));
    } else if (member.kind == SymbolKind::SpecifyBlock) {
      // `specify ... endspecify`: path delays / timing checks only,
      // no functional effect on the DUT's logic.
      logger.log(1,
                 "SystemVerilogEncoder: ignoring specify block (timing-only)");
    } else if (member.kind == SymbolKind::DefParam) {
      // `defparam inst.PARAM = value;`: a legacy parameter-override
      // mechanism with no functional-logic representation of its own
      // (the target instance is encoded with its own declared
      // defaults instead). The `#(...)` instantiation-time override
      // styles are fully supported; only this legacy spelling isn't.
      logger.log(1,
                 "SystemVerilogEncoder: ignoring defparam (parameter "
                 "overrides via `#(...)` at instantiation are supported; "
                 "`defparam` is not)");
    }
  });

  // Sequential and assertion-bearing blocks come second.
  walk_members(body, prefix_, [&](const Symbol & member) {
    if (member.kind == SymbolKind::ProceduralBlock) {
      auto & proc = member.as<ProceduralBlockSymbol>();
      switch (proc.procedureKind) {
        case ProceduralBlockKind::AlwaysFF: process_always_ff(proc); break;
        case ProceduralBlockKind::Initial: {
          // `initial forever @(...) body` is a legacy structural
          // spelling of `always @(...) body` -- redirect to the same
          // NEXT_STATE processing an always_ff block gets instead of
          // treating it as an initial-state constraint.
          if (auto * forever_body = as_forever_event_body(proc.getBody())) {
            process_next_state_body(*forever_body);
          } else {
            process_initial(proc);
          }
          break;
        }
        case ProceduralBlockKind::Always: {
          std::unordered_set<const Symbol *> targets;
          collect_nonblocking_targets(proc.getBody(), targets);
          if (!targets.empty()) {
            process_always_ff(proc);
          }
          break;
        }
        case ProceduralBlockKind::AlwaysLatch:
          // A level-sensitive latch's writes (blocking `=`, implicit
          // hold when a path doesn't reassign) are encoded exactly
          // like a register's (nonblocking `<=`, defaulting to itself
          // when not written) -- NEXT_STATE processing doesn't care
          // which assignment operator was used, only that writes
          // should become assign_next() targets.
          process_next_state_body(proc.getBody());
          break;
        default:
          // AlwaysComb handled above; Final, etc. skipped.
          break;
      }
    }
  });

  current_scope_ = saved_scope;
}

void SystemVerilogEncoder::process_always_ff(
    const slang::ast::ProceduralBlockSymbol & proc)
{
  process_next_state_body(proc.getBody());
}

void SystemVerilogEncoder::process_next_state_body(
    const slang::ast::Statement & body)
{
  symbol_table_.pending_next_updates().clear();

  // Use a null condition to represent "unconditional".
  Term true_term = solver_->make_term(true);
  const slang::ast::Expression * default_disable_expr =
      current_scope_ ? compilation_->getDefaultDisable(*current_scope_)
                     : nullptr;
  try {
    statement_encoder_.process_statement(
        body,
        StatementEncoder::StmtContext::NEXT_STATE,
        true_term,
        prefix_,
        default_disable_expr);
  }
  catch (const LoopControlSignal &) {
    // A compile-time-constant break/continue/disable is absorbed by a
    // matching enclosing ForLoop/named Block. A runtime-dependent one
    // is rejected right at the statement itself (it never becomes a
    // LoopControlSignal at all -- see StatementKind::Break/Continue/
    // Disable in statement.cpp). Only "no matching enclosing
    // construct anywhere" reaches this catch.
    throw PonoException(
        "SystemVerilogEncoder: break/continue/disable is only supported "
        "when its condition is a compile-time constant (e.g. depends only "
        "on already-unrolled for-loop counters)");
  }

  // Commit all pending next-state updates.
  for (auto & [state_term, next_expr] : symbol_table_.pending_next_updates()) {
    fts_.assign_next(state_term, next_expr);
    logger.log(2,
               "SystemVerilogEncoder: assign_next {} := ...",
               fts_.get_name(state_term));
  }
}

void SystemVerilogEncoder::process_always_comb(
    const slang::ast::ProceduralBlockSymbol & proc)
{
  symbol_table_.pending_comb_updates().clear();
  symbol_table_.pending_comb_aliased().clear();
  Term true_term = solver_->make_term(true);
  const slang::ast::Expression * default_disable_expr =
      current_scope_ ? compilation_->getDefaultDisable(*current_scope_)
                     : nullptr;
  try {
    statement_encoder_.process_statement(
        proc.getBody(),
        StatementEncoder::StmtContext::COMBINATIONAL,
        true_term,
        prefix_,
        default_disable_expr);
  }
  catch (const LoopControlSignal &) {
    throw PonoException(
        "SystemVerilogEncoder: break/continue/disable is only supported "
        "when its condition is a compile-time constant (e.g. depends only "
        "on already-unrolled for-loop counters)");
  }

  // Commit accumulated wire definitions via macro substitution.
  // Aliased entries belong in the parent's scope; everything else
  // uses the current prefix.
  for (auto & [sym, term] : symbol_table_.pending_comb_updates()) {
    string name;
    if (symbol_table_.pending_comb_aliased().count(sym)) {
      name = parent_prefix_.empty() ? string(sym->name)
                                    : parent_prefix_ + "." + string(sym->name);
    } else {
      name = symbol_table_.make_name(prefix_, string(sym->name));
    }
    symbol_table_.symbol_to_term()[sym] = term;
    fts_.name_term(name, term);
    logger.log(2, "SystemVerilogEncoder: always_comb (wire) {} := ...", name);
  }
  symbol_table_.pending_comb_updates().clear();
  symbol_table_.pending_comb_aliased().clear();
}

void SystemVerilogEncoder::process_initial(
    const slang::ast::ProceduralBlockSymbol & proc)
{
  Term true_term = solver_->make_term(true);
  const slang::ast::Expression * default_disable_expr =
      current_scope_ ? compilation_->getDefaultDisable(*current_scope_)
                     : nullptr;
  try {
    statement_encoder_.process_statement(proc.getBody(),
                                         StatementEncoder::StmtContext::INITIAL,
                                         true_term,
                                         prefix_,
                                         default_disable_expr);
  }
  catch (const LoopControlSignal &) {
    throw PonoException(
        "SystemVerilogEncoder: break/continue/disable is only supported "
        "when its condition is a compile-time constant (e.g. depends only "
        "on already-unrolled for-loop counters)");
  }
}

void SystemVerilogEncoder::process_always_comb_once(
    const slang::ast::ProceduralBlockSymbol & proc)
{
  if (!symbol_table_.processed_drivers().insert(&proc).second) return;
  process_always_comb(proc);
}

void SystemVerilogEncoder::process_continuous_assign_once(
    const slang::ast::ContinuousAssignSymbol & ca)
{
  if (!symbol_table_.processed_drivers().insert(&ca).second) return;
  process_continuous_assign(ca);
}

void SystemVerilogEncoder::process_continuous_assign(
    const slang::ast::ContinuousAssignSymbol & ca)
{
  using namespace slang::ast;

  auto & assign_expr = ca.getAssignment();
  if (assign_expr.kind != ExpressionKind::Assignment) {
    return;
  }

  auto & assign = assign_expr.as<AssignmentExpression>();
  auto & lhs_expr = assign.left();
  auto & rhs_expr = assign.right();

  // Concatenation-target LHS (`assign {hi, lo} = ...;`): unlike a
  // range-/element-select, this has more than one base symbol, so it
  // can't be represented as a single LValueDesc. Split the RHS across
  // each operand MSB-first (leftmost operand = most significant) and
  // recurse into process_continuous_assign_operand() once per operand,
  // exactly mirroring how a concatenation-target *port connection* is
  // already split into one OutputAliasSegment per operand.
  if (lhs_expr.kind == ExpressionKind::Concatenation) {
    // A concatenation's own type -- and so each slice written through
    // it -- is always unsigned per the LRM, regardless of the RHS
    // expression's own signedness: this is positional bit-splicing,
    // not a numeric value being widened.
    uint64_t total_w = lhs_expr.type->getBitWidth();
    if (total_w == 0) return;
    Term rhs_full = resize_to(
        solver_, expr_encoder_.expr_to_term(rhs_expr, prefix_), total_w, false);
    uint64_t covered = 0;
    for (auto * operand : lhs_expr.as<ConcatenationExpression>().operands()) {
      uint64_t seg_w = operand->type->getBitWidth();
      if (seg_w == 0 || covered + seg_w > total_w) break;
      uint64_t seg_hi = total_w - 1 - covered;
      uint64_t seg_lo = seg_hi - (seg_w - 1);
      process_continuous_assign_operand(
          *operand, slice_bits(solver_, rhs_full, seg_lo, seg_hi), false);
      covered += seg_w;
    }
    return;
  }

  process_continuous_assign_operand(
      lhs_expr,
      expr_encoder_.expr_to_term(rhs_expr, prefix_),
      rhs_expr.type->isSigned());
}

void SystemVerilogEncoder::process_continuous_assign_operand(
    const slang::ast::Expression & lhs_expr,
    const smt::Term & rhs_arg,
    bool rhs_signed)
{
  using namespace slang::ast;

  // Unlike the procedural-assignment path (which falls back to
  // process_dynamic_element_assign() for a genuinely dynamic-index
  // ElementSelect), a continuous assign has no such fallback -- a
  // nullopt here always means resolve_lvalue() couldn't statically
  // resolve this lvalue at all, so throw rather than silently
  // dropping the write.
  auto desc = resolve_lvalue(lhs_expr, expr_encoder_.eval_ctx());
  if (!desc) {
    throw PonoException(
        "SystemVerilogEncoder: unsupported continuous-assign lvalue "
        "(non-constant index?)");
  }
  const Symbol * base_sym = desc->base;
  bool aliased = symbol_table_.port_output_aliases().count(base_sym) > 0;

  Term rhs_full =
      resize_to(solver_, rhs_arg, desc->hi - desc->lo + 1, rhs_signed);

  // A concatenation-target output-port connection splits this one
  // write across several pieces, each with its own target
  // symbol/bit-range and its own slice of rhs_full; the common,
  // non-aliased (or singly-aliased) case is exactly one piece
  // spanning the whole write.
  auto pieces =
      symbol_table_.resolve_output_alias_pieces(base_sym, desc->lo, desc->hi);
  for (auto & piece : pieces) {
    const Symbol * sym = piece.sym;
    uint64_t lo = piece.target_lo;
    uint64_t hi = piece.target_hi;
    Term rhs = slice_bits(solver_, rhs_full, piece.rhs_lo, piece.rhs_hi);

    // Wire LHS: macro-substitute the *full-width* defining expression.
    // For a partial LHS (`assign arr[i] = ...`, or one element of an
    // instance array wired to a slice of a bus) we splice the slice
    // into whatever was previously stored under `sym`, creating a
    // fresh placeholder to splice into on the very first such write.
    if (symbol_table_.wire_symbols().count(sym)) {
      // Check against the symbol's own declared width, not the width
      // of whatever (possibly still-partial) term is already stored
      // under it -- otherwise a first write that starts at bit 0 but
      // doesn't cover the whole symbol gets mistaken for a full write,
      // corrupting the width of later, non-adjacent slice writes.
      uint64_t sym_w = sym->as<ValueSymbol>().getType().getBitWidth();
      bool full_write = (lo == 0 && hi + 1 == sym_w);
      Term new_term;
      if (full_write) {
        new_term = rhs;
      } else {
        new_term = replace_bits(
            solver_, symbol_table_.wire_seed_term(sym, prefix_), rhs, lo, hi);
      }
      symbol_table_.symbol_to_term()[sym] = new_term;
      // Only register the debug name on a full write: a wire spliced
      // together from several partial writes (e.g. separate sibling
      // instances each driving a different slice of one shared bus)
      // has no single write that owns its name, and re-naming it on
      // every partial write would collide with name_term()'s "one name,
      // one term" invariant as soon as two of those partial terms
      // differ.
      if (full_write) {
        string name;
        if (aliased) {
          name = parent_prefix_.empty()
                     ? string(sym->name)
                     : parent_prefix_ + "." + string(sym->name);
        } else {
          name = symbol_table_.make_name(prefix_, string(sym->name));
        }
        fts_.name_term(name, new_term);
        logger.log(2,
                   "SystemVerilogEncoder: continuous assign (wire) {} := ...",
                   name);
      }
      continue;
    }

    // Fallback: existing variable (e.g., output port reg, or a
    // partially-driven base that wasn't classified as a wire).
    // Constrain the appropriate slice via add_constraint (which
    // tolerates input vars in the term).
    auto it = symbol_table_.symbol_to_term().find(sym);
    if (it != symbol_table_.symbol_to_term().end()) {
      Term lhs_term = it->second;
      uint64_t base_w = lhs_term->get_sort()->get_width();
      bool full_write = (lo == 0 && hi == base_w - 1);
      Term lhs_slice = full_write
                           ? lhs_term
                           : solver_->make_term(Op(Extract, hi, lo), lhs_term);
      Term eq = solver_->make_term(Equal, lhs_slice, rhs);
      fts_.add_constraint(eq);
      logger.log(2,
                 "SystemVerilogEncoder: continuous assign {} = ...",
                 fts_.get_name(lhs_term));
    }
  }
}

void SystemVerilogEncoder::process_instance(
    const slang::ast::InstanceSymbol & inst)
{
  using namespace slang::ast;

  // `program ... endprogram`, instantiated like a module: a
  // verification-only construct (a testbench entry point, not
  // synthesizable DUT logic), so this instance and everything inside
  // it is simulation-only and has no functional-logic counterpart to
  // encode. (`interface` instances share SymbolKind::Instance with
  // ordinary modules too, but are a supported signal-bundle feature,
  // not a simulation-only one -- only DefinitionKind::Program is
  // skipped here.)
  if (inst.body.getDefinition().definitionKind == DefinitionKind::Program) {
    logger.log(1,
               "SystemVerilogEncoder: ignoring program instance '{}' "
               "(simulation-only construct)",
               string(inst.name));
    return;
  }

  // Switch the naming context so any state vars / wires declared
  // inside the child get hierarchical names like "<parent>.<inst>.<x>".
  // Also track the *parent's* prefix so wires redirected via
  // port_output_aliases_ (which live in the parent's scope) get
  // named correctly.
  string saved_prefix = prefix_;
  string saved_parent_prefix = parent_prefix_;
  parent_prefix_ = prefix_;
  // An instance-array element's name is empty (slang names only the
  // array itself); walk_members() has already pushed a "[i]"-suffixed
  // prefix for it, so don't append another separator here.
  prefix_ =
      inst.name.empty() ? saved_prefix : saved_prefix + "." + string(inst.name);

  // Bind the child's port-internal symbols to their parent-side
  // counterparts.  Inputs become parent-side terms (so reads inside
  // the child resolve via lookup_symbol).  Outputs become aliases (so
  // writes inside the child redirect to parent-side wires).  Save the
  // additions so we can undo them at the end of this call -- slang
  // may share an InstanceBody across multiple instantiations.
  std::vector<const Symbol *> output_aliases_added;
  std::vector<const Symbol *> input_terms_added;
  for (auto * pc : inst.getPortConnections()) {
    if (!pc) continue;
    if (pc->port.kind != SymbolKind::Port) continue;
    auto & port = pc->port.as<PortSymbol>();
    auto * conn_expr = pc->getExpression();
    if (!conn_expr) continue;
    auto * internal = port.internalSymbol;
    if (!internal) continue;

    bool is_output = (port.direction == ArgumentDirection::Out
                      || port.direction == ArgumentDirection::InOut);
    if (is_output) {
      // Output-port connections are wrapped in an Assignment whose
      // left-hand side is the parent-side expression.  For an
      // instance-array element, slang has already resolved this to
      // the correct constant-index slice of the parent-side bus
      // signal (e.g. `fifo_data_out[i]`), which resolve_lvalue()
      // decomposes into a base symbol and bit range just like any
      // other constant-index select.
      if (conn_expr->kind == ExpressionKind::Assignment) {
        conn_expr = &conn_expr->as<AssignmentExpression>().left();
      }
      uint64_t port_w = port.getType().getBitWidth();
      if (conn_expr->kind == ExpressionKind::Concatenation) {
        // `.port({hi, lo})`: split the port's bits across each
        // operand, MSB-first (leftmost operand = most significant),
        // one segment per operand.
        std::vector<OutputAliasSegment> segments;
        bool ok = port_w > 0;
        uint64_t covered = 0;
        for (auto * operand :
             conn_expr->as<ConcatenationExpression>().operands()) {
          auto odesc = resolve_lvalue(*operand, expr_encoder_.eval_ctx());
          if (!odesc) {
            ok = false;
            break;
          }
          uint64_t seg_w = odesc->hi - odesc->lo + 1;
          if (seg_w == 0 || covered + seg_w > port_w) {
            ok = false;
            break;
          }
          uint64_t seg_hi = port_w - 1 - covered;
          uint64_t seg_lo = seg_hi - (seg_w - 1);
          segments.push_back(
              { seg_lo, seg_hi, odesc->base, odesc->lo, odesc->hi });
          covered += seg_w;
        }
        if (ok && covered == port_w && !segments.empty()) {
          symbol_table_.port_output_aliases()[internal] = std::move(segments);
          output_aliases_added.push_back(internal);
        }
      } else {
        auto desc = resolve_lvalue(*conn_expr, expr_encoder_.eval_ctx());
        if (desc) {
          symbol_table_.port_output_aliases()[internal] = {
            { 0, port_w - 1, desc->base, desc->lo, desc->hi }
          };
          output_aliases_added.push_back(internal);
        }
      }
    } else {
      Term term = expr_encoder_.expr_to_term(*conn_expr, prefix_);
      term = resize_to(solver_,
                       term,
                       port.getType().getBitWidth(),
                       conn_expr->type->isSigned());
      symbol_table_.symbol_to_term()[internal] = term;
      input_terms_added.push_back(internal);
    }
  }

  // Pre-scan the child's blocking-assigned wires so they get classified
  // before declare_variables_internal runs. NB-assigned registers are
  // *not* pre-scanned here -- process_module()'s pre_scan_state_vars()
  // already classified every always_ff/always block in the whole
  // design tree, including this instance's, before any instance's
  // variables were declared.
  walk_members(inst.body, prefix_, [&](const Symbol & m) {
    if (m.kind != SymbolKind::ProceduralBlock) return;
    auto & proc = m.as<ProceduralBlockSymbol>();
    if (proc.procedureKind == ProceduralBlockKind::AlwaysComb) {
      symbol_table_.pre_scan_always_comb(
          proc.getBody(), proc, prefix_, parent_prefix_);
    } else if (proc.procedureKind == ProceduralBlockKind::Always) {
      std::unordered_set<const Symbol *> nb_targets;
      collect_nonblocking_targets(proc.getBody(), nb_targets);
      if (nb_targets.empty()) {
        symbol_table_.pre_scan_always_comb(
            proc.getBody(), proc, prefix_, parent_prefix_);
      }
    }
  });

  // Declare the child's internal (non-port) variables with the new
  // hierarchical prefix; ports are already bound through the
  // connection map above.
  declarer_.declare_variables_internal(inst.body, prefix_);

  // Combinational pass over child's body (and any sub-instances).
  walk_members(inst.body, prefix_, [&](const Symbol & m) {
    if (m.kind == SymbolKind::ContinuousAssign) {
      process_continuous_assign_once(m.as<ContinuousAssignSymbol>());
    } else if (m.kind == SymbolKind::ProceduralBlock) {
      auto & proc = m.as<ProceduralBlockSymbol>();
      if (proc.procedureKind == ProceduralBlockKind::AlwaysComb) {
        process_always_comb_once(proc);
      } else if (proc.procedureKind == ProceduralBlockKind::Always) {
        std::unordered_set<const Symbol *> targets;
        collect_nonblocking_targets(proc.getBody(), targets);
        if (targets.empty()) {
          process_always_comb_once(proc);
        }
      }
    } else if (m.kind == SymbolKind::Instance) {
      process_instance(m.as<InstanceSymbol>());
    }
  });

  // Sequential / initial pass.
  walk_members(inst.body, prefix_, [&](const Symbol & m) {
    if (m.kind != SymbolKind::ProceduralBlock) return;
    auto & proc = m.as<ProceduralBlockSymbol>();
    switch (proc.procedureKind) {
      case ProceduralBlockKind::AlwaysFF: process_always_ff(proc); break;
      case ProceduralBlockKind::Initial: {
        if (auto * forever_body = as_forever_event_body(proc.getBody())) {
          process_next_state_body(*forever_body);
        } else {
          process_initial(proc);
        }
        break;
      }
      case ProceduralBlockKind::Always: {
        std::unordered_set<const Symbol *> targets;
        collect_nonblocking_targets(proc.getBody(), targets);
        if (!targets.empty()) process_always_ff(proc);
        break;
      }
      case ProceduralBlockKind::AlwaysLatch:
        process_next_state_body(proc.getBody());
        break;
      default: break;
    }
  });

  // Restore context: undo the per-instance bindings so that a sibling
  // (or repeated) instantiation of the same module can be processed
  // cleanly.
  for (auto * sym : output_aliases_added) {
    symbol_table_.port_output_aliases().erase(sym);
  }
  for (auto * sym : input_terms_added) {
    symbol_table_.symbol_to_term().erase(sym);
  }
  prefix_ = saved_prefix;
  parent_prefix_ = saved_parent_prefix;
}

}  // namespace pono
