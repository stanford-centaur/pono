/*!
 * \file instance_encoder.cpp
 * \brief Per-instance pass: continuous assigns, procedural blocks, instances.
 * \author Áron Ricardo Perez-Lopez
 * \date 2026
 * \copyright See the LICENSE file in the top-level source directory.
 */
#include "frontends/systemverilog/instance_encoder.h"

#include <unordered_set>
#include <vector>

#include "frontends/systemverilog/ast_helpers.h"
#include "frontends/systemverilog/bit_utils.h"
#include "frontends/systemverilog/declarer.h"
#include "frontends/systemverilog/expr_encoder.h"
#include "frontends/systemverilog/statement_encoder.h"
#include "slang/ast/Compilation.h"
#include "slang/ast/EvalContext.h"
#include "slang/ast/Expression.h"
#include "slang/ast/SemanticFacts.h"
#include "slang/ast/Statement.h"
#include "slang/ast/Symbol.h"
#include "slang/ast/expressions/AssignmentExpressions.h"
#include "slang/ast/expressions/MiscExpressions.h"
#include "slang/ast/expressions/OperatorExpressions.h"
#include "slang/ast/expressions/SelectExpressions.h"
#include "slang/ast/symbols/BlockSymbols.h"
#include "slang/ast/symbols/CheckerSymbols.h"
#include "slang/ast/symbols/CompilationUnitSymbols.h"
#include "slang/ast/symbols/InstanceSymbols.h"
#include "slang/ast/symbols/MemberSymbols.h"
#include "slang/ast/symbols/PortSymbols.h"
#include "slang/ast/types/AllTypes.h"
#include "slang/ast/types/Type.h"
#include "smt-switch/smt.h"
#include "utils/exceptions.h"
#include "utils/logger.h"

using namespace smt;
using namespace std;

namespace pono {

InstanceEncoder::InstanceEncoder(SymbolTable & symbol_table,
                                 Declarer & declarer,
                                 StatementEncoder & statement_encoder,
                                 ExprEncoder & expr_encoder,
                                 FunctionalTransitionSystem & fts,
                                 const smt::SmtSolver & solver)
    : symbol_table_(symbol_table),
      declarer_(declarer),
      statement_encoder_(statement_encoder),
      expr_encoder_(expr_encoder),
      fts_(fts),
      solver_(solver)
{
  symbol_table_.set_driver_resolver(*this);
}

void InstanceEncoder::bind_compilation(slang::ast::Compilation & compilation)
{
  compilation_ = &compilation;
}

namespace {

// Re-express alias segments given in generic-stream coordinates as
// segments of the `<<` re-ordered port. Both coordinate systems run
// low bit to high within a segment, and the re-ordering moves whole
// blocks without turning them over, so each run of bits that stays
// contiguous in both is one output segment -- there are just more of
// them than operands now, one per piece of an operand that a block
// boundary cuts.
//
// Walking bit by bit and coalescing, rather than intersecting
// intervals, keeps this honest about the short block: it sits at the
// opposite end from the one it was cut from, and is the case an
// interval calculation gets subtly wrong.
std::vector<OutputAliasSegment> reblock_stream_segments(
    const std::vector<OutputAliasSegment> & generic,
    uint64_t port_w,
    uint64_t slice)
{
  std::vector<OutputAliasSegment> out;
  for (const OutputAliasSegment & seg : generic) {
    bool open = false;
    uint64_t run_port_lo = 0, run_port_hi = 0;
    uint64_t run_target_lo = 0, run_target_hi = 0;
    for (uint64_t g = seg.port_lo; g <= seg.port_hi; ++g) {
      uint64_t p = stream_reorder_bit(g, port_w, slice);
      uint64_t t = seg.target_lo + (g - seg.port_lo);
      if (open && p == run_port_hi + 1 && t == run_target_hi + 1) {
        run_port_hi = p;
        run_target_hi = t;
        continue;
      }
      if (open) {
        out.push_back({ run_port_lo,
                        run_port_hi,
                        seg.target,
                        run_target_lo,
                        run_target_hi });
      }
      run_port_lo = run_port_hi = p;
      run_target_lo = run_target_hi = t;
      open = true;
    }
    if (open) {
      out.push_back({ run_port_lo,
                      run_port_hi,
                      seg.target,
                      run_target_lo,
                      run_target_hi });
    }
  }
  return out;
}

}  // namespace

void InstanceEncoder::process_assignments(const slang::ast::Scope & body,
                                          const string & prefix,
                                          const string & parent_prefix)
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

  // walk_members() takes the prefix by mutable reference (updating it
  // while descending into generate-for/instance-array child scopes,
  // then restoring it), so it needs its own local copy rather than
  // `prefix` itself.
  string walk_prefix = prefix;

  // Combinational definitions are processed first so that wires have a
  // term assigned in symbol_to_term_ before any always_ff or initial
  // block (or assertion) tries to reference them.  Child instances are
  // walked here too -- a child's continuous assigns / always_comb
  // blocks may drive parent-side wires that downstream parent code
  // references.
  walk_members(body, walk_prefix, [&](const Symbol & member) {
    if (member.kind == SymbolKind::ContinuousAssign) {
      process_continuous_assign_once(
          member.as<ContinuousAssignSymbol>(), walk_prefix, parent_prefix);
    } else if (member.kind == SymbolKind::ProceduralBlock) {
      auto & proc = member.as<ProceduralBlockSymbol>();
      if (proc.procedureKind == ProceduralBlockKind::AlwaysComb) {
        process_always_comb_once(proc, walk_prefix, parent_prefix);
      } else if (proc.procedureKind == ProceduralBlockKind::Always) {
        std::unordered_set<const Symbol *> targets;
        collect_nonblocking_targets(proc.getBody(), targets);
        // An edge-sensitive `always` is a register block even when it
        // writes only with `=`, so it belongs to the sequential walk
        // below; without the edge test a blocking-only one matches
        // neither walk and is dropped.
        if (targets.empty() && !is_edge_triggered(proc.getBody())) {
          process_always_comb_once(proc, walk_prefix, parent_prefix);
        }
      }
    } else if (member.kind == SymbolKind::Instance) {
      process_instance(member.as<InstanceSymbol>(), walk_prefix, parent_prefix);
    } else if (member.kind == SymbolKind::CheckerInstance) {
      process_checker_instance(
          member.as<CheckerInstanceSymbol>(), walk_prefix, parent_prefix);
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
  walk_prefix = prefix;
  walk_members(body, walk_prefix, [&](const Symbol & member) {
    if (member.kind == SymbolKind::ProceduralBlock) {
      auto & proc = member.as<ProceduralBlockSymbol>();
      switch (proc.procedureKind) {
        case ProceduralBlockKind::AlwaysFF:
          process_always_ff(proc, walk_prefix);
          break;
        case ProceduralBlockKind::Initial: {
          // `initial forever @(...) body` is a legacy structural
          // spelling of `always @(...) body` -- redirect to the same
          // NEXT_STATE processing an always_ff block gets instead of
          // treating it as an initial-state constraint.
          if (auto * forever_body = as_forever_event_body(proc.getBody())) {
            process_next_state_body(*forever_body, walk_prefix);
          } else {
            process_initial(proc, walk_prefix);
          }
          break;
        }
        case ProceduralBlockKind::Always: {
          std::unordered_set<const Symbol *> targets;
          collect_nonblocking_targets(proc.getBody(), targets);
          if (!targets.empty() || is_edge_triggered(proc.getBody())) {
            process_always_ff(proc, walk_prefix);
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
          process_next_state_body(proc.getBody(), walk_prefix);
          break;
        default:
          // AlwaysComb handled above; Final, etc. skipped.
          break;
      }
    }
  });

  current_scope_ = saved_scope;
}

void InstanceEncoder::process_always_ff(
    const slang::ast::ProceduralBlockSymbol & proc, const string & prefix)
{
  process_next_state_body(proc.getBody(), prefix);
}

void InstanceEncoder::process_next_state_body(
    const slang::ast::Statement & body, const string & prefix)
{
  symbol_table_.pending_next_updates().clear();
  symbol_table_.blocking_next_written().clear();

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
        prefix,
        default_disable_expr);
  }
  catch (const LoopControlSignal &) {
    // A compile-time-constant break/continue/disable is absorbed by a
    // matching enclosing ForLoop/named Block. A runtime-dependent one
    // is rejected right at the statement itself (it never becomes a
    // LoopControlSignal at all -- see StatementKind::Break/Continue/
    // Disable in statement_encoder.cpp). Only "no matching enclosing
    // construct anywhere" reaches this catch.
    throw PonoException(
        "SystemVerilogEncoder: break/continue/disable is only supported "
        "when its condition is a compile-time constant (e.g. depends only "
        "on already-unrolled for-loop counters)");
  }

  // A local that got storage but that no path in this block writes
  // holds its value for good; without this it would have no
  // next-state function and so be free every cycle, which is a looser
  // reading than "unknown, but the same unknown".
  std::unordered_set<const slang::ast::Symbol *> holds;
  collect_hold_locals(body, holds);
  for (auto * sym : holds) {
    auto sit = symbol_table_.symbol_to_term().find(sym);
    if (sit == symbol_table_.symbol_to_term().end()) continue;
    if (symbol_table_.pending_next_updates().count(sit->second)) continue;
    symbol_table_.pending_next_updates()[sit->second] = sit->second;
  }

  // Commit all pending next-state updates.
  for (auto & [state_term, next_expr] : symbol_table_.pending_next_updates()) {
    fts_.assign_next(state_term, next_expr);
    logger.log(2,
               "SystemVerilogEncoder: assign_next {} := ...",
               fts_.get_name(state_term));
  }
}

bool InstanceEncoder::holds_previous_value(const Term & target,
                                           const Term & value)
{
  Sort sort = target->get_sort();
  string base = "__latch_probe_" + std::to_string(latch_probe_counter_++);
  Term a = solver_->make_symbol(base + "_a", sort);
  Term b = solver_->make_symbol(base + "_b", sort);
  UnorderedTermMap with_a{ { target, a } };
  UnorderedTermMap with_b{ { target, b } };
  Term differ = solver_->make_term(Distinct,
                                   solver_->substitute(value, with_a),
                                   solver_->substitute(value, with_b));
  bool can_differ = true;
  try {
    solver_->push();
    solver_->assert_formula(differ);
    can_differ = !solver_->check_sat().is_unsat();
    solver_->pop();
  }
  catch (const std::exception &) {
    return true;
  }
  return can_differ;
}

void InstanceEncoder::process_always_comb(
    const slang::ast::ProceduralBlockSymbol & proc,
    const string & prefix,
    const string & parent_prefix)
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
        prefix,
        default_disable_expr);
  }
  catch (const LoopControlSignal &) {
    throw PonoException(
        "SystemVerilogEncoder: break/continue/disable is only supported "
        "when its condition is a compile-time constant (e.g. depends only "
        "on already-unrolled for-loop counters)");
  }

  // Commit each accumulated definition. A wire is macro-substituted
  // (aliased entries belong in the parent's scope, everything else
  // uses the current prefix); a non-wire target already has a term of
  // its own, so it takes a single constraint equating that term to
  // the value the whole block computed -- one per symbol, since a
  // constraint per write would bind the same term several times over.
  for (auto & [sym, term] : symbol_table_.pending_comb_updates()) {
    if (!symbol_table_.wire_symbols().count(sym)) {
      auto sit = symbol_table_.symbol_to_term().find(sym);
      if (sit == symbol_table_.symbol_to_term().end()) {
        throw PonoException("SystemVerilogEncoder: always_comb writes '"
                            + string(sym->name)
                            + "', which has no declared term");
      }
      // A path that writes nothing leaves the accumulated value
      // falling back to the symbol's own term, and the scan that
      // spotted that is syntactic: it counts a `case` with no
      // `default` as leaving a path open even when the arms between
      // them cover every value. Believing it there would turn plain
      // combinational logic into a register, delaying it a cycle --
      // which proves properties the design does not have. So ask
      // whether that fallback can actually be reached.
      if (symbol_table_.latch_symbols().count(sym)
          && holds_previous_value(sit->second, term)) {
        // It can, so the target really does keep its old value on
        // those paths -- which a same-cycle equality cannot say, and
        // a next-state update says exactly. Warned about for the same
        // reason a synthesis tool warns: in an `always_comb` a latch
        // is rarely intended.
        logger.log(0,
                   "SystemVerilogEncoder: '{}' is assigned on only some "
                   "paths through a combinational block, so it holds its "
                   "value on the rest -- modeled as the latch this infers",
                   string(sym->name));
        fts_.assign_next(sit->second, term);
        continue;
      }
      fts_.add_constraint(solver_->make_term(Equal, sit->second, term));
      logger.log(2,
                 "SystemVerilogEncoder: always_comb (reg) {} := ...",
                 fts_.get_name(sit->second));
      continue;
    }
    string name;
    if (symbol_table_.pending_comb_aliased().count(sym)) {
      name = parent_prefix.empty() ? string(sym->name)
                                   : parent_prefix + "." + string(sym->name);
    } else {
      name = symbol_table_.make_name(prefix, string(sym->name));
    }
    symbol_table_.symbol_to_term()[sym] = term;
    fts_.name_term(name, term);
    logger.log(2, "SystemVerilogEncoder: always_comb (wire) {} := ...", name);
  }
  symbol_table_.pending_comb_updates().clear();
  symbol_table_.pending_comb_aliased().clear();
}

void InstanceEncoder::process_initial(
    const slang::ast::ProceduralBlockSymbol & proc, const string & prefix)
{
  symbol_table_.pending_comb_updates().clear();
  Term true_term = solver_->make_term(true);
  const slang::ast::Expression * default_disable_expr =
      current_scope_ ? compilation_->getDefaultDisable(*current_scope_)
                     : nullptr;
  try {
    statement_encoder_.process_statement(proc.getBody(),
                                         StatementEncoder::StmtContext::INITIAL,
                                         true_term,
                                         prefix,
                                         default_disable_expr);
  }
  catch (const LoopControlSignal &) {
    throw PonoException(
        "SystemVerilogEncoder: break/continue/disable is only supported "
        "when its condition is a compile-time constant (e.g. depends only "
        "on already-unrolled for-loop counters)");
  }

  // One constraint per symbol, pinning it to the value the whole
  // block left it with -- a constraint per write would bind the same
  // term several times over and contradict itself.
  for (auto & [sym, term] : symbol_table_.pending_comb_updates()) {
    auto sit = symbol_table_.symbol_to_term().find(sym);
    if (sit == symbol_table_.symbol_to_term().end()) {
      throw PonoException("SystemVerilogEncoder: initial block writes '"
                          + string(sym->name)
                          + "', which has no declared term");
    }
    Term init_eq = solver_->make_term(Equal, sit->second, term);
    if (!fts_.only_curr(init_eq)) {
      // An initial block runs before any input has a meaning, so a
      // value it depends on has to be part of the design's state.
      // The core rejects anything else, in terms that say nothing
      // about the block this came from.
      throw PonoException(
          "SystemVerilogEncoder: the initial value of '" + string(sym->name)
          + "' depends on an input, which has no value at time 0");
    }
    fts_.constrain_init(init_eq);
    logger.log(2,
               "SystemVerilogEncoder: initial {} := ...",
               fts_.get_name(sit->second));
  }
  symbol_table_.pending_comb_updates().clear();
}

void InstanceEncoder::process_always_comb_once(
    const slang::ast::ProceduralBlockSymbol & proc,
    const string & prefix,
    const string & parent_prefix)
{
  if (!symbol_table_.processed_drivers().insert(&proc).second) return;
  process_always_comb(proc, prefix, parent_prefix);
}

void InstanceEncoder::process_continuous_assign_once(
    const slang::ast::ContinuousAssignSymbol & ca,
    const string & prefix,
    const string & parent_prefix)
{
  if (!symbol_table_.processed_drivers().insert(&ca).second) return;
  process_continuous_assign(ca, prefix, parent_prefix);
}

void InstanceEncoder::resolve_continuous_assign(
    const slang::ast::ContinuousAssignSymbol & ca,
    const string & prefix,
    const string & parent_prefix)
{
  process_continuous_assign_once(ca, prefix, parent_prefix);
}

void InstanceEncoder::resolve_always_comb(
    const slang::ast::ProceduralBlockSymbol & proc,
    const string & prefix,
    const string & parent_prefix)
{
  process_always_comb_once(proc, prefix, parent_prefix);
}

void InstanceEncoder::process_continuous_assign(
    const slang::ast::ContinuousAssignSymbol & ca,
    const string & prefix,
    const string & parent_prefix)
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
  // A streaming concatenation target (`assign {>>{hi, lo}} = ...;`)
  // splits the same way; what it adds is that the source is consumed
  // from its most significant end and then un-re-ordered.
  if (lhs_expr.kind == ExpressionKind::Concatenation
      || lhs_expr.kind == ExpressionKind::Streaming) {
    bool streaming = lhs_expr.kind == ExpressionKind::Streaming;
    std::vector<const Expression *> operands;
    uint64_t total_w = 0;
    if (streaming) {
      auto & sc = lhs_expr.as<StreamingConcatenationExpression>();
      total_w = sc.getBitstreamWidth();
      for (auto & stream : sc.streams()) {
        if (stream.withExpr) {
          throw PonoException(
              "SystemVerilogEncoder: a `with` range in a streaming "
              "concatenation target is not supported");
        }
        operands.push_back(stream.operand);
      }
    } else {
      total_w = lhs_expr.type->getBitWidth();
      for (auto * operand : lhs_expr.as<ConcatenationExpression>().operands()) {
        operands.push_back(operand);
      }
    }
    if (total_w == 0) return;
    // A concatenation's own type -- and so each slice written through
    // it -- is always unsigned per the LRM, regardless of the RHS
    // expression's own signedness: this is positional bit-splicing,
    // not a numeric value being widened.
    Term rhs_full = expr_encoder_.expr_to_term(rhs_expr, prefix);
    if (streaming) {
      uint64_t rhs_w = rhs_full->get_sort()->get_width();
      if (rhs_w < total_w) {
        throw PonoException(
            "SystemVerilogEncoder: a streaming-concatenation target needs "
            + std::to_string(total_w) + " bits but the source supplies only "
            + std::to_string(rhs_w));
      }
      rhs_full = slice_bits(solver_, rhs_full, rhs_w - total_w, rhs_w - 1);
      rhs_full = stream_unreorder(
          solver_,
          rhs_full,
          lhs_expr.as<StreamingConcatenationExpression>().getSliceSize());
    } else {
      rhs_full = resize_to(solver_, rhs_full, total_w, false);
    }
    uint64_t covered = 0;
    for (const Expression * operand : operands) {
      uint64_t seg_w = value_width(*operand->type);
      if (seg_w == 0 || covered + seg_w > total_w) break;
      uint64_t seg_hi = total_w - 1 - covered;
      uint64_t seg_lo = seg_hi - (seg_w - 1);
      process_continuous_assign_operand(
          *operand,
          slice_bits(solver_, rhs_full, seg_lo, seg_hi),
          false,
          prefix,
          parent_prefix);
      covered += seg_w;
    }
    return;
  }

  process_continuous_assign_operand(
      lhs_expr,
      expr_encoder_.expr_to_term(rhs_expr, prefix),
      rhs_expr.type->isSigned(),
      prefix,
      parent_prefix);
}

void InstanceEncoder::process_continuous_assign_operand(
    const slang::ast::Expression & lhs_expr,
    const smt::Term & rhs_arg,
    bool rhs_signed,
    const string & prefix,
    const string & parent_prefix)
{
  using namespace slang::ast;

  // An unpacked array is not a bit range, so resolve_lvalue() has
  // nothing to say about it. Driving one continuously is a
  // constraint rather than an assignment: the whole array equals the
  // right-hand side, or one element of it does. Elements left
  // undriven stay free, which is what an undriven net is.
  {
    const Expression * target = &lhs_expr;
    const Expression * index = nullptr;
    if (target->kind == ExpressionKind::ElementSelect) {
      auto & sel = target->as<ElementSelectExpression>();
      if (sel.value().type->getCanonicalType().kind
          == SymbolKind::FixedSizeUnpackedArrayType) {
        index = &sel.selector();
        target = &sel.value();
      }
    }
    bool whole_array = !index
                       && target->type->getCanonicalType().kind
                              == SymbolKind::FixedSizeUnpackedArrayType;
    if (whole_array && target->kind == ExpressionKind::NamedValue) {
      // A wire takes its value from its driver, so for a whole-array
      // one this assign *is* the definition -- there is no term yet
      // to constrain.
      const Symbol * sym =
          &canonicalize_signal_alias(target->as<NamedValueExpression>().symbol);
      if (!symbol_table_.symbol_to_term().count(sym)) {
        symbol_table_.symbol_to_term()[sym] = rhs_arg;
        return;
      }
    }
    if (index || whole_array) {
      Term array = expr_encoder_.expr_to_term(*target, prefix);
      Term driven = array;
      if (index) {
        auto & arr =
            target->type->getCanonicalType().as<FixedSizeUnpackedArrayType>();
        UnpackedArrayInfo info = unpacked_array_info(solver_, arr);
        Term idx = resize_to(solver_,
                             expr_encoder_.expr_to_term(*index, prefix),
                             info.index_width,
                             index->type->isSigned());
        driven = solver_->make_term(
            Select, array, normalize_array_index(solver_, idx, info));
      }
      fts_.add_constraint(solver_->make_term(Equal, driven, rhs_arg));
      return;
    }
  }

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
      uint64_t sym_w = value_width(sym->as<ValueSymbol>().getType());
      bool full_write = (lo == 0 && hi + 1 == sym_w);
      Term new_term;
      if (full_write) {
        new_term = rhs;
      } else {
        new_term = replace_bits(
            solver_, symbol_table_.wire_seed_term(sym, prefix), rhs, lo, hi);
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
          name = parent_prefix.empty()
                     ? string(sym->name)
                     : parent_prefix + "." + string(sym->name);
        } else {
          name = symbol_table_.make_name(prefix, string(sym->name));
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
    if (it == symbol_table_.symbol_to_term().end()) {
      // A hierarchical continuous-assign target with no declared term
      // yet means this assign's source position precedes the child
      // instance's own declaration in the same scope (declaration is
      // interleaved with, and ordered by, source position) -- or, for
      // a plain (non-hierarchical) target, that it reaches into a
      // child instance's internal (non-port) signal from outside that
      // instance's own scope at all. Neither is real synthesizable
      // RTL (module ports are the only sanctioned cross-instance
      // wiring mechanism), so throw rather than silently dropping the
      // write and leaving the target fully unconstrained.
      throw PonoException(
          "SystemVerilogEncoder: unsupported continuous-assign target '"
          + string(sym->name)
          + "' (hierarchical reference into a child instance's internal "
            "signal, or forward reference to an instance declared later "
            "in the same scope?)");
    }
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

void InstanceEncoder::process_instance(const slang::ast::InstanceSymbol & inst,
                                       const string & prefix,
                                       const string & parent_prefix,
                                       bool assertions_only)
{
  using namespace slang::ast;

  // `program ... endprogram`, instantiated like a module: a testbench
  // entry point rather than synthesizable DUT logic. Its stimulus is
  // deliberately left unencoded -- a program drives the DUT along one
  // particular scenario, and pinning the inputs to it would leave
  // every other input sequence unexplored while still reporting a
  // proof. An assertion written inside one is not stimulus, though:
  // it is an ordinary property, and dropping it means reporting that
  // proof over fewer properties than were written. So encode those
  // and nothing else. (`interface` instances share
  // SymbolKind::Instance with ordinary modules too, but are a
  // supported signal-bundle feature, not a simulation-only one --
  // only DefinitionKind::Program is treated this way.)
  if (inst.body.getDefinition().definitionKind == DefinitionKind::Program) {
    logger.log(1,
               "SystemVerilogEncoder: encoding only the concurrent "
               "assertions of program instance '{}'; its stimulus is "
               "simulation-only",
               string(inst.name));
    assertions_only = true;
  }

  // Compute the child's own hierarchical prefix, and track the
  // *parent's* prefix (this call's own `prefix`) so wires redirected
  // via port_output_aliases_ (which live in the parent's scope) get
  // named correctly. Plain local variables, not a mutated-and-restored
  // shared member: each recursive process_instance() call down the
  // instance tree gets its own copy on the call stack.
  //
  // An instance-array element's name is empty (slang names only the
  // array itself); walk_members() has already pushed a "[i]"-suffixed
  // prefix for it, so don't append another separator here.
  string child_prefix =
      inst.name.empty() ? prefix : prefix + "." + string(inst.name);
  const string & child_parent_prefix = prefix;

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
    // An explicit port names its internal signal through an
    // expression rather than the internalSymbol field; reading that
    // field directly skipped the connection without a word, leaving
    // the child's side of it an unconstrained variable.
    auto * internal = port_internal_symbol(port);
    if (!internal) continue;

    bool is_output = (port.direction == ArgumentDirection::Out
                      || port.direction == ArgumentDirection::InOut);
    if (port.direction == ArgumentDirection::InOut) {
      // Only the child's drive is modelled. That is exactly right
      // while nothing else drives the net, which covers an `inout`
      // used in one direction at a time; what it cannot represent is
      // the parent driving back, since resolving two drivers onto
      // one net is not something this encoder has any notion of.
      logger.log(0,
                 "SystemVerilogEncoder: inout port '{}' of instance '{}' "
                 "is modeled as an output -- the child drives it and "
                 "anything driving it from outside is not modeled",
                 string(port.name),
                 string(inst.name));
    }
    if (is_output && assertions_only) {
      // An output of a program is stimulus, which is not being
      // encoded. Registering an alias redirects the parent-side
      // target to writes that will now never arrive, leaving it
      // undeclared; leaving the port alone keeps it an ordinary
      // unconstrained signal instead.
      continue;
    }
    if (is_output) {
      // Output-port connections are wrapped in an Assignment whose
      // left-hand side is the parent-side expression.  For an
      // instance-array element, slang has already resolved this to
      // the correct constant-index slice of the parent-side bus
      // signal (e.g. `fifo_data_out[i]`), which resolve_lvalue()
      // decomposes into a base symbol and bit range just like any
      // other constant-index select.
      //
      // Unlike the procedural-assignment path (which falls back to
      // process_dynamic_element_assign() for a genuinely dynamic-index
      // ElementSelect), a port connection has no such fallback -- it's
      // a structural, elaboration-time binding, not a per-cycle write,
      // so there's no mux to build. A resolve_lvalue() failure here
      // always means the connection expression is unsupported (e.g. a
      // non-constant index), so throw rather than silently dropping
      // the port's write and leaving the target fully unconstrained.
      if (conn_expr->kind == ExpressionKind::Assignment) {
        conn_expr = &conn_expr->as<AssignmentExpression>().left();
      }
      uint64_t port_w = value_width(port.getType());
      if (conn_expr->kind == ExpressionKind::Concatenation
          || conn_expr->kind == ExpressionKind::Streaming) {
        // `.port({hi, lo})`: split the port's bits across each
        // operand, MSB-first (leftmost operand = most significant),
        // one segment per operand.
        std::vector<const Expression *> operands;
        // Zero for `>>` and for a plain concatenation, both of which
        // re-order nothing.
        uint64_t slice = 0;
        if (conn_expr->kind == ExpressionKind::Streaming) {
          auto & sc = conn_expr->as<StreamingConcatenationExpression>();
          slice = sc.getSliceSize();
          for (auto & stream : sc.streams()) {
            if (stream.withExpr) {
              throw PonoException(
                  "SystemVerilogEncoder: a `with` range in a streaming "
                  "concatenation connected to output/inout port '"
                  + string(port.name) + "' is not supported");
            }
            operands.push_back(stream.operand);
          }
        } else {
          for (auto * operand :
               conn_expr->as<ConcatenationExpression>().operands()) {
            operands.push_back(operand);
          }
        }
        std::vector<OutputAliasSegment> segments;
        uint64_t covered = 0;
        for (const Expression * operand : operands) {
          auto odesc = resolve_lvalue(*operand, expr_encoder_.eval_ctx());
          if (!odesc) {
            throw PonoException(
                "SystemVerilogEncoder: unsupported concatenation-target "
                "output/inout port connection for port '"
                + string(port.name)
                + "' (non-constant index in a concatenation operand)");
          }
          uint64_t seg_w = odesc->hi - odesc->lo + 1;
          if (seg_w == 0 || covered + seg_w > port_w) {
            throw PonoException(
                "SystemVerilogEncoder: output/inout port connection for "
                "port '" + string(port.name)
                + "' has a concatenation width that doesn't fit the "
                "port width");
          }
          uint64_t seg_hi = port_w - 1 - covered;
          uint64_t seg_lo = seg_hi - (seg_w - 1);
          segments.push_back(
              { seg_lo, seg_hi, odesc->base, odesc->lo, odesc->hi });
          covered += seg_w;
        }
        if (covered != port_w || segments.empty()) {
          throw PonoException(
              "SystemVerilogEncoder: output/inout port connection for "
              "port '" + string(port.name)
              + "' concatenation total width does not match the port "
              "width");
        }
        // The operands were laid out against the generic stream. A
        // `<<` re-orders that into the port's own bits, so a segment
        // spanning a block boundary is no longer one port range and
        // has to be cut where the blocks are.
        if (slice != 0) {
          segments = reblock_stream_segments(segments, port_w, slice);
        }
        symbol_table_.port_output_aliases()[internal] = std::move(segments);
        output_aliases_added.push_back(internal);
      } else if (port.getType().getCanonicalType().kind
                 == SymbolKind::FixedSizeUnpackedArrayType) {
        // A whole array has no bits to splice across alias segments,
        // so the child's port variable and the parent's array are
        // simply the same term, and the child's writes land on it
        // directly. The parent's array is a state variable because
        // the child drives it; nothing in the parent does, so no
        // other pass would have made it one.
        const Symbol * target = find_lhs_base(*conn_expr);
        if (!target) {
          throw PonoException(
              "SystemVerilogEncoder: the array connected to output port '"
              + string(port.name) + "' is not a plain variable");
        }
        Sort port_sort = type_to_sort(solver_, port.getType());
        auto * target_value = target->as_if<ValueSymbol>();
        if (!target_value) {
          throw PonoException(
              "SystemVerilogEncoder: the array connected to output port '"
              + string(port.name) + "' is not a plain variable");
        }
        Sort target_sort = type_to_sort(solver_, target_value->getType());

        // A slice cannot share the parent's term -- the two have
        // different lengths. Give the port its own array and tie it
        // to the parent's, element by element: the count is known at
        // elaboration, so this is a fixed handful of equalities
        // rather than anything dynamic.
        if (target_sort != port_sort) {
          if (conn_expr->kind != ExpressionKind::RangeSelect) {
            throw PonoException(
                "SystemVerilogEncoder: output port '" + string(port.name)
                + "' is connected to part of '" + string(target->name)
                + "' that is not a contiguous element range");
          }
          auto & rs = conn_expr->as<RangeSelectExpression>();
          auto & base_ct = rs.value().type->getCanonicalType();
          auto lc = rs.left().eval(expr_encoder_.eval_ctx());
          auto rc = rs.right().eval(expr_encoder_.eval_ctx());
          if (rs.getSelectionKind() != RangeSelectionKind::Simple
              || base_ct.kind != SymbolKind::FixedSizeUnpackedArrayType
              || !lc.isInteger() || !rc.isInteger()) {
            throw PonoException(
                "SystemVerilogEncoder: output port '" + string(port.name)
                + "' is connected to a slice of '" + string(target->name)
                + "' whose bounds are not elaboration-time constants");
          }
          auto lv = lc.integer().as<int64_t>();
          auto rv = rc.integer().as<int64_t>();
          if (!lv || !rv) {
            throw PonoException("SystemVerilogEncoder: output port '"
                                + string(port.name)
                                + "' has a slice bound that does not fit");
          }
          auto & base_arr = base_ct.as<FixedSizeUnpackedArrayType>();
          UnpackedArrayInfo base_info = unpacked_array_info(solver_, base_arr);
          UnpackedArrayInfo port_info =
              unpacked_array_info(solver_,
                                  port.getType()
                                      .getCanonicalType()
                                      .as<FixedSizeUnpackedArrayType>());
          // Normalized indices count up from the declared range's
          // lower bound, so the slice's own low end is the offset
          // into the parent however either range is written.
          int64_t low = std::min(*lv, *rv);
          int64_t offset = low - base_arr.range.lower();
          if (offset < 0
              || static_cast<uint64_t>(offset) + port_info.depth
                     > base_info.depth) {
            throw PonoException(
                "SystemVerilogEncoder: output port '" + string(port.name)
                + "' is connected to a slice of '" + string(target->name)
                + "' that runs outside it");
          }

          auto pit = symbol_table_.symbol_to_term().find(target);
          Term parent_term;
          if (pit != symbol_table_.symbol_to_term().end()) {
            parent_term = pit->second;
          } else {
            parent_term = fts_.make_statevar(
                symbol_table_.make_name(prefix, string(target->name)),
                target_sort);
            symbol_table_.symbol_to_term()[target] = parent_term;
            symbol_table_.state_var_symbols().insert(target);
            symbol_table_.wire_symbols().erase(target);
          }
          Term child_term = fts_.make_statevar(
              symbol_table_.make_name(child_prefix, string(internal->name)),
              port_sort);
          Sort port_idx = solver_->make_sort(BV, port_info.index_width);
          Sort base_idx = solver_->make_sort(BV, base_info.index_width);
          for (uint64_t k = 0; k < port_info.depth; ++k) {
            fts_.add_constraint(solver_->make_term(
                Equal,
                solver_->make_term(
                    Select, child_term, solver_->make_term(k, port_idx)),
                solver_->make_term(Select,
                                   parent_term,
                                   solver_->make_term(offset + k, base_idx))));
          }
          symbol_table_.symbol_to_term()[internal] = child_term;
          input_terms_added.push_back(internal);
          continue;
        }
        auto it = symbol_table_.symbol_to_term().find(target);
        Term shared;
        if (it != symbol_table_.symbol_to_term().end()) {
          shared = it->second;
          if (shared->get_sort() != port_sort) {
            throw PonoException(
                "SystemVerilogEncoder: the array connected to output port '"
                + string(port.name) + "' has sort "
                + shared->get_sort()->to_string() + ", but the port is "
                + port_sort->to_string());
          }
        } else {
          shared = fts_.make_statevar(
              symbol_table_.make_name(prefix, string(target->name)), port_sort);
          symbol_table_.symbol_to_term()[target] = shared;
          symbol_table_.state_var_symbols().insert(target);
          symbol_table_.wire_symbols().erase(target);
        }
        // Only the child's own symbol is undone afterwards; the
        // parent keeps its array.
        symbol_table_.symbol_to_term()[internal] = shared;
        input_terms_added.push_back(internal);
      } else {
        auto desc = resolve_lvalue(*conn_expr, expr_encoder_.eval_ctx());
        if (!desc) {
          throw PonoException(
              "SystemVerilogEncoder: unsupported output/inout port "
              "connection for port '"
              + string(port.name) + "' (non-constant index?)");
        }
        symbol_table_.port_output_aliases()[internal] = {
          { 0, port_w - 1, desc->base, desc->lo, desc->hi }
        };
        output_aliases_added.push_back(internal);
      }
    } else {
      Term term = expr_encoder_.expr_to_term(*conn_expr, prefix);
      if (term->get_sort()->get_sort_kind() == ARRAY) {
        // An unpacked array is passed whole -- there is no width to
        // reconcile, only the requirement that both sides agree.
        Sort port_sort = type_to_sort(solver_, port.getType());
        if (term->get_sort() != port_sort) {
          throw PonoException(
              "SystemVerilogEncoder: the array connected to port '"
              + string(port.name) + "' has sort "
              + term->get_sort()->to_string() + ", but the port is "
              + port_sort->to_string());
        }
      } else {
        term = resize_to(solver_,
                         term,
                         port.getType().getBitWidth(),
                         conn_expr->type->isSigned());
      }
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
  string walk_prefix = child_prefix;
  walk_members(inst.body, walk_prefix, [&](const Symbol & m) {
    if (m.kind != SymbolKind::ProceduralBlock) return;
    auto & proc = m.as<ProceduralBlockSymbol>();
    if (proc.procedureKind == ProceduralBlockKind::AlwaysComb) {
      symbol_table_.pre_scan_always_comb(
          proc.getBody(), proc, walk_prefix, child_parent_prefix);
    } else if (proc.procedureKind == ProceduralBlockKind::Always) {
      std::unordered_set<const Symbol *> nb_targets;
      collect_nonblocking_targets(proc.getBody(), nb_targets);
      if (nb_targets.empty()) {
        symbol_table_.pre_scan_always_comb(
            proc.getBody(), proc, walk_prefix, child_parent_prefix);
      }
    }
  });

  // Declare the child's internal (non-port) variables with the new
  // hierarchical prefix; ports are already bound through the
  // connection map above.
  declarer_.declare_variables_internal(inst.body, child_prefix);

  // Combinational pass over child's body (and any sub-instances).
  // A module-scope concurrent assertion reaches the encoder as an
  // `Always` block carrying nothing but the assertion, so this is the
  // pass that picks one up -- which is why `assertions_only` filters
  // here rather than skipping the walk.
  walk_prefix = child_prefix;
  walk_members(inst.body, walk_prefix, [&](const Symbol & m) {
    if (assertions_only) {
      if (m.kind != SymbolKind::ProceduralBlock) return;
      auto & proc = m.as<ProceduralBlockSymbol>();
      if (is_concurrent_assertion_only(proc.getBody())) {
        process_always_comb_once(proc, walk_prefix, child_parent_prefix);
      }
      return;
    }
    if (m.kind == SymbolKind::ContinuousAssign) {
      process_continuous_assign_once(
          m.as<ContinuousAssignSymbol>(), walk_prefix, child_parent_prefix);
    } else if (m.kind == SymbolKind::ProceduralBlock) {
      auto & proc = m.as<ProceduralBlockSymbol>();
      if (proc.procedureKind == ProceduralBlockKind::AlwaysComb) {
        process_always_comb_once(proc, walk_prefix, child_parent_prefix);
      } else if (proc.procedureKind == ProceduralBlockKind::Always) {
        std::unordered_set<const Symbol *> targets;
        collect_nonblocking_targets(proc.getBody(), targets);
        if (targets.empty()) {
          process_always_comb_once(proc, walk_prefix, child_parent_prefix);
        }
      }
    } else if (m.kind == SymbolKind::Instance) {
      process_instance(
          m.as<InstanceSymbol>(), walk_prefix, child_parent_prefix);
    } else if (m.kind == SymbolKind::CheckerInstance) {
      process_checker_instance(
          m.as<CheckerInstanceSymbol>(), walk_prefix, child_parent_prefix);
    } else if (m.kind == SymbolKind::SpecifyBlock) {
      logger.log(1,
                 "SystemVerilogEncoder: ignoring specify block (timing-only)");
    } else if (m.kind == SymbolKind::DefParam) {
      logger.log(1,
                 "SystemVerilogEncoder: ignoring defparam (parameter "
                 "overrides via `#(...)` at instantiation are supported; "
                 "`defparam` is not)");
    }
  });

  // Sequential / initial pass. Everything it would encode is
  // stimulus, so there is nothing here for a program instance.
  if (assertions_only) return;
  walk_prefix = child_prefix;
  walk_members(inst.body, walk_prefix, [&](const Symbol & m) {
    if (m.kind != SymbolKind::ProceduralBlock) return;
    auto & proc = m.as<ProceduralBlockSymbol>();
    switch (proc.procedureKind) {
      case ProceduralBlockKind::AlwaysFF:
        process_always_ff(proc, walk_prefix);
        break;
      case ProceduralBlockKind::Initial: {
        if (auto * forever_body = as_forever_event_body(proc.getBody())) {
          process_next_state_body(*forever_body, walk_prefix);
        } else {
          process_initial(proc, walk_prefix);
        }
        break;
      }
      case ProceduralBlockKind::Always: {
        std::unordered_set<const Symbol *> targets;
        collect_nonblocking_targets(proc.getBody(), targets);
        if (!targets.empty()) process_always_ff(proc, walk_prefix);
        break;
      }
      case ProceduralBlockKind::AlwaysLatch:
        process_next_state_body(proc.getBody(), walk_prefix);
        break;
      default: break;
    }
  });

  // Undo the per-instance bindings so that a sibling (or repeated)
  // instantiation of the same module can be processed cleanly.
  for (auto * sym : output_aliases_added) {
    symbol_table_.port_output_aliases().erase(sym);
  }
  for (auto * sym : input_terms_added) {
    symbol_table_.symbol_to_term().erase(sym);
  }
}

void InstanceEncoder::process_checker_instance(
    const slang::ast::CheckerInstanceSymbol & ci,
    const string & prefix,
    const string & parent_prefix)
{
  using namespace slang::ast;

  string child_prefix =
      ci.name.empty() ? prefix : prefix + "." + string(ci.name);
  const string & child_parent_prefix = prefix;

  // Unlike a module instance, a checker's formal (`AssertionPortSymbol`)
  // ports are resolved by slang itself at elaboration time: a
  // reference to a formal inside the checker's body already binds
  // directly to the actual argument's own symbol (e.g. the caller's
  // `clk` net), not to a distinct checker-local copy. So there is no
  // port-binding step to do here at all, unlike process_instance()'s
  // alias/input-term setup above.
  //
  // A checker's own genuinely local state -- a `Variable`/`Net`
  // declared directly in its body, or a combinational blocking-
  // assigned target -- is otherwise exactly like a module instance's:
  // pre-scan its blocking-assigned wires (nonblocking-target state
  // vars were already classified up front, whole-tree, by
  // SystemVerilogEncoder::process_module()'s pre_scan_state_vars()
  // call, which now also recurses into checker instances), then
  // declare its internal variables under the checker instance's own
  // hierarchical prefix, mirroring process_instance() above exactly.
  string walk_prefix = child_prefix;
  walk_members(ci.body, walk_prefix, [&](const Symbol & m) {
    if (m.kind != SymbolKind::ProceduralBlock) return;
    auto & proc = m.as<ProceduralBlockSymbol>();
    if (proc.procedureKind == ProceduralBlockKind::AlwaysComb) {
      symbol_table_.pre_scan_always_comb(
          proc.getBody(), proc, walk_prefix, child_parent_prefix);
    } else if (proc.procedureKind == ProceduralBlockKind::Always) {
      std::unordered_set<const Symbol *> nb_targets;
      collect_nonblocking_targets(proc.getBody(), nb_targets);
      if (nb_targets.empty()) {
        symbol_table_.pre_scan_always_comb(
            proc.getBody(), proc, walk_prefix, child_parent_prefix);
      }
    }
  });

  declarer_.declare_variables_internal(ci.body, child_prefix);

  process_assignments(ci.body, child_prefix, child_parent_prefix);
}

}  // namespace pono
