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
 **/

#include "frontends/systemverilog_encoder.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "core/fts.h"
#include "slang/ast/Compilation.h"
#include "slang/ast/EvalContext.h"
#include "slang/ast/Expression.h"
#include "slang/ast/Scope.h"
#include "slang/ast/SemanticFacts.h"
#include "slang/ast/Statement.h"
#include "slang/ast/Symbol.h"
#include "slang/ast/expressions/AssertionExpr.h"
#include "slang/ast/expressions/AssignmentExpressions.h"
#include "slang/ast/expressions/ConversionExpression.h"
#include "slang/ast/expressions/LiteralExpressions.h"
#include "slang/ast/expressions/MiscExpressions.h"
#include "slang/ast/expressions/Operator.h"
#include "slang/ast/expressions/OperatorExpressions.h"
#include "slang/ast/expressions/SelectExpressions.h"
#include "slang/ast/statements/ConditionalStatements.h"
#include "slang/ast/statements/LoopStatements.h"
#include "slang/ast/statements/MiscStatements.h"
#include "slang/ast/symbols/BlockSymbols.h"
#include "slang/ast/symbols/CompilationUnitSymbols.h"
#include "slang/ast/symbols/InstanceSymbols.h"
#include "slang/ast/symbols/MemberSymbols.h"
#include "slang/ast/symbols/ParameterSymbols.h"
#include "slang/ast/symbols/PortSymbols.h"
#include "slang/ast/symbols/VariableSymbols.h"
#include "slang/ast/types/Type.h"
#include "slang/diagnostics/DiagnosticEngine.h"
#include "slang/numeric/SVInt.h"
#include "slang/syntax/SyntaxTree.h"
#include "smt-switch/smt.h"
#include "utils/exceptions.h"
#include "utils/logger.h"

using namespace smt;
using namespace std;

namespace pono {

// Forward declarations for helpers defined in an anonymous namespace
// later in this file.
namespace {
void collect_nonblocking_targets(
    const slang::ast::Statement & stmt,
    std::unordered_set<const slang::ast::Symbol *> & targets);
// Collects bases of blocking assignments inside `stmt`.  `full` gets
// only NamedValue / HierarchicalValue (full-width) LHS bases;
// `partial` gets bases reached through ElementSelect / RangeSelect.
// Symbols that appear in both sets are partially written somewhere
// and so must be classified as state vars to keep add_invar slice
// constraints valid.
void collect_blocking_targets(
    const slang::ast::Statement & stmt,
    std::unordered_set<const slang::ast::Symbol *> & full,
    std::unordered_set<const slang::ast::Symbol *> & partial);

// Identifies the base ValueSymbol underlying a (possibly nested)
// bit/range-select LHS.  Returns nullptr if the LHS shape isn't
// supported by the encoder (e.g. concatenation LHS).
const slang::ast::Symbol * find_lhs_base(const slang::ast::Expression & lhs);

// Describes an LHS slice: which base symbol gets written and at what
// bit range.  For a NamedValue (or HierarchicalValue) LHS this is the
// full range [0, base_w-1]; for nested ElementSelects of constant
// indices the range narrows accordingly while base_w stays the full
// base bit width.
struct LValueDesc
{
  const slang::ast::Symbol * base;
  uint64_t lo;
  uint64_t hi;
  uint64_t base_w;
};

std::optional<LValueDesc> resolve_lvalue(const slang::ast::Expression & lhs,
                                         slang::ast::EvalContext & ctx);
}  // namespace

// ============================================================================
// Member iteration helpers
// ============================================================================

template <typename Fn>
void SystemVerilogEncoder::walk_members(const slang::ast::Scope & scope,
                                        Fn && fn)
{
  using namespace slang::ast;
  for (auto & m : scope.members()) {
    if (m.kind == SymbolKind::GenerateBlockArray) {
      // Generate-for: walk each instantiated entry, pushing a
      // bracket-indexed prefix so per-iteration variables get
      // unique hierarchical names like "<top>.ctr[0].count".
      auto & arr = m.as<GenerateBlockArraySymbol>();
      std::string saved_prefix = prefix_;
      std::string arr_name = std::string(arr.name);
      if (arr_name.empty()) arr_name = arr.getExternalName();
      for (auto * entry : arr.entries) {
        if (!entry || entry->isUninstantiated) continue;
        std::string idx_str;
        if (entry->arrayIndex) {
          auto idx = *entry->arrayIndex;
          idx.setSigned(false);
          idx_str =
              idx.toString(slang::LiteralBase::Decimal, /*includeBase=*/false);
        } else {
          idx_str = std::to_string(entry->constructIndex);
        }
        prefix_ = saved_prefix + "." + arr_name + "[" + idx_str + "]";
        walk_members(*entry, fn);
      }
      prefix_ = saved_prefix;
    } else if (m.kind == SymbolKind::GenerateBlock) {
      // Generate-if / generate-case: a single block scope.  Push
      // its name (or slang's synthesized "genblkN") as the suffix.
      auto & gb = m.as<GenerateBlockSymbol>();
      if (gb.isUninstantiated) continue;
      std::string saved_prefix = prefix_;
      std::string block_name = std::string(gb.name);
      if (block_name.empty()) block_name = gb.getExternalName();
      prefix_ = saved_prefix + "." + block_name;
      walk_members(gb, fn);
      prefix_ = saved_prefix;
    } else {
      fn(m);
    }
  }
}

slang::ast::EvalContext & SystemVerilogEncoder::eval_ctx()
{
  if (!eval_ctx_) {
    // Use the slang compilation root as the AST context scope; we
    // only use the eval context's locals stack to bind procedural
    // loop counters, so the lookup location doesn't matter.
    slang::ast::ASTContext ast_ctx(compilation_->getRoot(),
                                   slang::ast::LookupLocation::min);
    eval_ctx_ = std::make_unique<slang::ast::EvalContext>(ast_ctx);
    eval_ctx_->pushEmptyFrame();
  }
  return *eval_ctx_;
}

// ============================================================================
// Construction / destruction
// ============================================================================

SystemVerilogEncoder::SystemVerilogEncoder(string filename,
                                           FunctionalTransitionSystem & fts)
    : fts_(fts), solver_(fts.solver())
{
  encode(filename);
}

SystemVerilogEncoder::~SystemVerilogEncoder() = default;

// ============================================================================
// Top-level encoding pipeline
// ============================================================================

void SystemVerilogEncoder::encode(const string & filename)
{
  // Parse the SystemVerilog source file.
  // SyntaxTree::fromFile returns an expected<shared_ptr<SyntaxTree>, ...>.
  auto tree_result = slang::syntax::SyntaxTree::fromFile(filename);
  if (!tree_result) {
    throw PonoException("SystemVerilogEncoder: failed to parse file: "
                        + filename);
  }
  auto tree = std::move(tree_result).value();

  // Create a compilation and elaborate the design.
  compilation_ = make_unique<slang::ast::Compilation>();
  compilation_->addSyntaxTree(tree);

  // Check for diagnostics (errors/warnings).
  // Force full elaboration before checking diagnostics.
  auto & diagnostics = compilation_->getAllDiagnostics();
  bool has_errors = false;
  for (size_t i = 0; i < diagnostics.size(); i++) {
    if (diagnostics[i].isError()) {
      has_errors = true;
      break;
    }
  }
  if (has_errors) {
    // Format diagnostics for the error message.
    auto * sm = compilation_->getSourceManager();
    string diag_messages;
    if (sm) {
      diag_messages = slang::DiagnosticEngine::reportAll(*sm, diagnostics);
    } else {
      diag_messages = "(unable to format diagnostics)";
    }
    throw PonoException("SystemVerilogEncoder: errors in SystemVerilog file:\n"
                        + diag_messages);
  }

  // Walk the top-level instances.
  auto & root = compilation_->getRoot();
  auto top_instances = root.topInstances;
  if (top_instances.empty()) {
    throw PonoException(
        "SystemVerilogEncoder: no top-level module instances found in "
        + filename);
  }

  // Process the first top-level module.
  // Multi-top designs could be supported by iterating, but for model
  // checking we typically have a single top module.
  process_module(*top_instances[0]);
}

// ============================================================================
// Module processing
// ============================================================================

void SystemVerilogEncoder::process_module(
    const slang::ast::InstanceSymbol & inst)
{
  string inst_name(inst.name);
  logger.log(1, "SystemVerilogEncoder: processing module {}", inst_name);
  prefix_ = inst_name;

  auto & body = inst.body;

  // First pass: identify state variable symbols by scanning always_ff blocks
  // for non-blocking assignment targets, before declaring variables.
  walk_members(body, [&](const slang::ast::Symbol & member) {
    if (member.kind == slang::ast::SymbolKind::ProceduralBlock) {
      auto & proc = member.as<slang::ast::ProceduralBlockSymbol>();
      if (proc.procedureKind == slang::ast::ProceduralBlockKind::AlwaysFF
          || proc.procedureKind == slang::ast::ProceduralBlockKind::Always) {
        // Pre-scan: find all non-blocking assignment targets.
        pre_scan_always_ff(proc.getBody());
      }
    }
  });

  // Second pre-pass: identify combinational wire symbols from
  // always_comb blocks, legacy `always` blocks without non-blocking
  // assignments, and continuous-assign LHS values.  Wires are not
  // independent variables -- they will be macro-substituted with their
  // defining expressions, so we must skip declaring them as input vars.
  walk_members(body, [&](const slang::ast::Symbol & member) {
    if (member.kind == slang::ast::SymbolKind::ProceduralBlock) {
      auto & proc = member.as<slang::ast::ProceduralBlockSymbol>();
      if (proc.procedureKind == slang::ast::ProceduralBlockKind::AlwaysComb) {
        pre_scan_always_comb(proc.getBody());
      } else if (proc.procedureKind
                 == slang::ast::ProceduralBlockKind::Always) {
        // Legacy always: combinational iff it has no non-blocking
        // assignments to identify it as sequential.
        std::unordered_set<const slang::ast::Symbol *> nb_targets;
        collect_nonblocking_targets(proc.getBody(), nb_targets);
        if (nb_targets.empty()) {
          pre_scan_always_comb(proc.getBody());
        }
      }
    } else if (member.kind == slang::ast::SymbolKind::ContinuousAssign) {
      auto & ca = member.as<slang::ast::ContinuousAssignSymbol>();
      auto & ae = ca.getAssignment();
      if (ae.kind == slang::ast::ExpressionKind::Assignment) {
        auto & lhs = ae.as<slang::ast::AssignmentExpression>().left();
        if (lhs.kind == slang::ast::ExpressionKind::NamedValue) {
          auto * sym = &lhs.as<slang::ast::NamedValueExpression>().symbol;
          if (!state_var_symbols_.count(sym)) {
            wire_symbols_.insert(sym);
          }
        } else if (auto * base = find_lhs_base(lhs)) {
          // Partial-LHS continuous assign (`assign arr[i] = ...`):
          // the base needs to be a state var so process_continuous_assign
          // can constrain the slice via add_invar.
          if (!wire_symbols_.count(base)) {
            state_var_symbols_.insert(base);
          }
        }
      }
    } else if (member.kind == slang::ast::SymbolKind::Instance) {
      pre_scan_instance(member.as<slang::ast::InstanceSymbol>());
    }
  });

  // Third pass: declare state vars, inputs, and output-port aliases.
  // Wire symbols are skipped here -- they get their defining term
  // assigned during combinational-assignment processing.
  declare_variables(body);

  // Fourth pass: process behavioral code and continuous assignments.
  process_assignments(body);
}

// ============================================================================
// Pre-scan: identify state variables from always_ff blocks
// ============================================================================

// We need a helper that's not declared in the header since it's
// implementation-only. We'll use a local recursive function via a
// namespace-scope helper.
namespace {

/// Helper to iterate over sub-statements of a BlockStatement body.
/// The body is a single Statement; if it is a StatementList, iterate its
/// children, otherwise visit the single statement directly.
template <typename Func>
void for_each_stmt_in_block(const slang::ast::BlockStatement & block,
                            Func && func)
{
  auto & body = block.body;
  if (body.kind == slang::ast::StatementKind::List) {
    auto & list = body.as<slang::ast::StatementList>();
    for (auto * s : list.list) {
      func(*s);
    }
  } else {
    func(body);
  }
}

// Collect blocking-assignment targets, separating full-width LHSes
// (wire candidates) from partial LHSes (which must be state vars so
// the assignment handler's add_invar slice constraints are valid).
void collect_blocking_targets(
    const slang::ast::Statement & stmt,
    std::unordered_set<const slang::ast::Symbol *> & full,
    std::unordered_set<const slang::ast::Symbol *> & partial)
{
  using namespace slang::ast;

  switch (stmt.kind) {
    case StatementKind::ExpressionStatement: {
      auto & es = stmt.as<ExpressionStatement>();
      auto & expr = es.expr;
      if (expr.kind == ExpressionKind::Assignment) {
        auto & assign = expr.as<AssignmentExpression>();
        auto & lhs = assign.left();
        if (lhs.kind == ExpressionKind::NamedValue) {
          full.insert(&lhs.as<NamedValueExpression>().symbol);
        } else if (lhs.kind == ExpressionKind::HierarchicalValue) {
          full.insert(&lhs.as<HierarchicalValueExpression>().symbol);
        } else if (auto * base = find_lhs_base(lhs)) {
          partial.insert(base);
        }
      }
      break;
    }
    case StatementKind::Block: {
      auto & block = stmt.as<BlockStatement>();
      auto & body = block.body;
      if (body.kind == StatementKind::List) {
        auto & list = body.as<StatementList>();
        for (auto * s : list.list) collect_blocking_targets(*s, full, partial);
      } else {
        collect_blocking_targets(body, full, partial);
      }
      break;
    }
    case StatementKind::Conditional: {
      auto & cond = stmt.as<ConditionalStatement>();
      collect_blocking_targets(cond.ifTrue, full, partial);
      if (cond.ifFalse) collect_blocking_targets(*cond.ifFalse, full, partial);
      break;
    }
    case StatementKind::Case: {
      auto & cs = stmt.as<CaseStatement>();
      for (auto & item : cs.items)
        collect_blocking_targets(*item.stmt, full, partial);
      if (cs.defaultCase)
        collect_blocking_targets(*cs.defaultCase, full, partial);
      break;
    }
    case StatementKind::Timed: {
      auto & ts = stmt.as<TimedStatement>();
      collect_blocking_targets(ts.stmt, full, partial);
      break;
    }
    case StatementKind::ForLoop: {
      // Recurse into the body so writes inside the
      // (compile-time-unrolled) loop are seen during pre-scan.
      auto & loop = stmt.as<ForLoopStatement>();
      collect_blocking_targets(loop.body, full, partial);
      break;
    }
    default: break;
  }
}

void collect_nonblocking_targets(
    const slang::ast::Statement & stmt,
    std::unordered_set<const slang::ast::Symbol *> & targets)
{
  using namespace slang::ast;

  switch (stmt.kind) {
    case StatementKind::ExpressionStatement: {
      auto & es = stmt.as<ExpressionStatement>();
      auto & expr = es.expr;
      if (expr.kind == ExpressionKind::Assignment) {
        auto & assign = expr.as<AssignmentExpression>();
        if (assign.isNonBlocking()) {
          // The LHS of a non-blocking assignment is a state variable;
          // for partial writes (`arr[i] <= ...`) we classify the base.
          if (auto * base = find_lhs_base(assign.left())) {
            targets.insert(base);
          }
        }
      }
      break;
    }
    case StatementKind::Block: {
      auto & block = stmt.as<BlockStatement>();
      for_each_stmt_in_block(block, [&](const Statement & s) {
        collect_nonblocking_targets(s, targets);
      });
      break;
    }
    case StatementKind::Conditional: {
      auto & cond = stmt.as<ConditionalStatement>();
      collect_nonblocking_targets(cond.ifTrue, targets);
      if (cond.ifFalse) {
        collect_nonblocking_targets(*cond.ifFalse, targets);
      }
      break;
    }
    case StatementKind::Case: {
      auto & cs = stmt.as<CaseStatement>();
      for (auto & item : cs.items) {
        collect_nonblocking_targets(*item.stmt, targets);
      }
      if (cs.defaultCase) {
        collect_nonblocking_targets(*cs.defaultCase, targets);
      }
      break;
    }
    case StatementKind::Timed: {
      auto & ts = stmt.as<TimedStatement>();
      collect_nonblocking_targets(ts.stmt, targets);
      break;
    }
    case StatementKind::ForLoop: {
      // Recurse into the body so NB-assigned registers inside the
      // (compile-time-unrolled) loop are seen during pre-scan.
      auto & loop = stmt.as<ForLoopStatement>();
      collect_nonblocking_targets(loop.body, targets);
      break;
    }
    default:
      // Other statement types: nothing to extract.
      break;
  }
}

const slang::ast::Symbol * find_lhs_base(const slang::ast::Expression & lhs)
{
  using namespace slang::ast;
  switch (lhs.kind) {
    case ExpressionKind::NamedValue:
      return &lhs.as<NamedValueExpression>().symbol;
    case ExpressionKind::HierarchicalValue:
      return &lhs.as<HierarchicalValueExpression>().symbol;
    case ExpressionKind::ElementSelect:
      return find_lhs_base(lhs.as<ElementSelectExpression>().value());
    case ExpressionKind::RangeSelect:
      return find_lhs_base(lhs.as<RangeSelectExpression>().value());
    default: return nullptr;
  }
}

std::optional<LValueDesc> resolve_lvalue(const slang::ast::Expression & lhs,
                                         slang::ast::EvalContext & ctx)
{
  using namespace slang::ast;
  switch (lhs.kind) {
    case ExpressionKind::NamedValue: {
      auto * sym = &lhs.as<NamedValueExpression>().symbol;
      uint64_t w = lhs.type->getBitWidth();
      if (w == 0) return std::nullopt;
      return LValueDesc{ sym, 0, w - 1, w };
    }
    case ExpressionKind::HierarchicalValue: {
      auto * sym = &lhs.as<HierarchicalValueExpression>().symbol;
      uint64_t w = lhs.type->getBitWidth();
      if (w == 0) return std::nullopt;
      return LValueDesc{ sym, 0, w - 1, w };
    }
    case ExpressionKind::ElementSelect: {
      auto & sel = lhs.as<ElementSelectExpression>();
      auto inner = resolve_lvalue(sel.value(), ctx);
      if (!inner) return std::nullopt;
      auto idx_cv = sel.selector().eval(ctx);
      if (!idx_cv.isInteger()) return std::nullopt;
      auto idx_opt = idx_cv.integer().as<uint64_t>();
      if (!idx_opt) return std::nullopt;
      uint64_t idx = *idx_opt;
      uint64_t elem_w = lhs.type->getBitWidth();
      if (elem_w == 0) return std::nullopt;
      uint64_t lo = inner->lo + idx * elem_w;
      uint64_t hi = lo + elem_w - 1;
      if (hi > inner->hi) return std::nullopt;
      return LValueDesc{ inner->base, lo, hi, inner->base_w };
    }
    default: return std::nullopt;
  }
}

}  // anonymous namespace

// This is a private helper called from process_module; it is not in the
// header to avoid exposing slang::ast::Statement there.
void SystemVerilogEncoder::pre_scan_always_ff(
    const slang::ast::Statement & body)
{
  collect_nonblocking_targets(body, state_var_symbols_);
}

void SystemVerilogEncoder::pre_scan_always_comb(
    const slang::ast::Statement & body)
{
  // Collect blocking-assign LHS targets.  Bases written full-width
  // become wires (macro-substituted); bases written through bit /
  // range selects become state vars instead so the assignment
  // handler can use slice-equality add_invar constraints (which
  // requires the term to be a state var, not an input var).  Mixed
  // full+partial writes also go to state vars to keep the splice
  // semantics correct.
  std::unordered_set<const slang::ast::Symbol *> full, partial;
  collect_blocking_targets(body, full, partial);
  for (auto * sym : full) {
    if (state_var_symbols_.count(sym)) continue;
    if (partial.count(sym)) {
      state_var_symbols_.insert(sym);
    } else {
      wire_symbols_.insert(sym);
    }
  }
  for (auto * sym : partial) {
    if (!wire_symbols_.count(sym)) {
      state_var_symbols_.insert(sym);
    }
  }
}

void SystemVerilogEncoder::pre_scan_instance(
    const slang::ast::InstanceSymbol & inst)
{
  using namespace slang::ast;

  // Each output (or inout) port of the child is driven by the child's
  // logic.  In the parent's view the connected expression must be a
  // simple NamedValue (e.g. `child c (.sum(y))`) -- the named symbol
  // becomes a wire.  More elaborate connections (slices,
  // concatenations, etc.) are not yet supported.
  for (auto * pc : inst.getPortConnections()) {
    if (!pc) continue;
    if (pc->port.kind != SymbolKind::Port) continue;
    auto & port = pc->port.as<PortSymbol>();
    if (port.direction != ArgumentDirection::Out
        && port.direction != ArgumentDirection::InOut)
      continue;
    auto * conn_expr = pc->getExpression();
    if (!conn_expr) continue;
    // Slang wraps an output-port connection as an Assignment whose
    // left-hand side is the parent-side expression being driven.
    // Unwrap it; the child's logic effectively writes through to the
    // LHS.
    if (conn_expr->kind == ExpressionKind::Assignment) {
      conn_expr = &conn_expr->as<AssignmentExpression>().left();
    }
    if (conn_expr->kind != ExpressionKind::NamedValue) continue;
    auto * parent_sym = &conn_expr->as<NamedValueExpression>().symbol;
    if (!state_var_symbols_.count(parent_sym)) {
      wire_symbols_.insert(parent_sym);
    }
  }

  // Recurse into nested instances so any wires further down the
  // hierarchy are visible to declare_variables.
  walk_members(inst.body, [&](const Symbol & m) {
    if (m.kind == SymbolKind::Instance) {
      pre_scan_instance(m.as<InstanceSymbol>());
    }
  });
}

// ============================================================================
// Variable declaration (first pass)
// ============================================================================

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
      // skip declaring a separate term for it.
      if (port_output_aliases_.count(&var)) return;

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

// ============================================================================
// Assignment processing (second pass)
// ============================================================================

void SystemVerilogEncoder::process_assignments(
    const slang::ast::InstanceBodySymbol & body)
{
  using namespace slang::ast;

  // Combinational definitions are processed first so that wires have a
  // term assigned in symbol_to_term_ before any always_ff or initial
  // block (or assertion) tries to reference them.  Child instances are
  // walked here too -- a child's continuous assigns / always_comb
  // blocks may drive parent-side wires that downstream parent code
  // references.
  walk_members(body, [&](const Symbol & member) {
    if (member.kind == SymbolKind::ContinuousAssign) {
      process_continuous_assign(member.as<ContinuousAssignSymbol>());
    } else if (member.kind == SymbolKind::ProceduralBlock) {
      auto & proc = member.as<ProceduralBlockSymbol>();
      if (proc.procedureKind == ProceduralBlockKind::AlwaysComb) {
        process_always_comb(proc);
      } else if (proc.procedureKind == ProceduralBlockKind::Always) {
        std::unordered_set<const Symbol *> targets;
        collect_nonblocking_targets(proc.getBody(), targets);
        if (targets.empty()) {
          process_always_comb(proc);
        }
      }
    } else if (member.kind == SymbolKind::Instance) {
      process_instance(member.as<InstanceSymbol>());
    }
  });

  // Sequential and assertion-bearing blocks come second.
  walk_members(body, [&](const Symbol & member) {
    if (member.kind == SymbolKind::ProceduralBlock) {
      auto & proc = member.as<ProceduralBlockSymbol>();
      switch (proc.procedureKind) {
        case ProceduralBlockKind::AlwaysFF: process_always_ff(proc); break;
        case ProceduralBlockKind::Initial: process_initial(proc); break;
        case ProceduralBlockKind::Always: {
          std::unordered_set<const Symbol *> targets;
          collect_nonblocking_targets(proc.getBody(), targets);
          if (!targets.empty()) {
            process_always_ff(proc);
          }
          break;
        }
        default:
          // AlwaysComb handled above; AlwaysLatch, Final, etc. skipped.
          break;
      }
    }
  });
}

void SystemVerilogEncoder::process_always_ff(
    const slang::ast::ProceduralBlockSymbol & proc)
{
  pending_next_updates_.clear();

  // Use a null condition to represent "unconditional".
  Term true_term = solver_->make_term(true);
  process_statement(proc.getBody(), StmtContext::NEXT_STATE, true_term);

  // Commit all pending next-state updates.
  for (auto & [state_term, next_expr] : pending_next_updates_) {
    fts_.assign_next(state_term, next_expr);
    logger.log(2,
               "SystemVerilogEncoder: assign_next {} := ...",
               fts_.get_name(state_term));
  }
}

void SystemVerilogEncoder::process_always_comb(
    const slang::ast::ProceduralBlockSymbol & proc)
{
  pending_comb_updates_.clear();
  pending_comb_aliased_.clear();
  Term true_term = solver_->make_term(true);
  process_statement(proc.getBody(), StmtContext::COMBINATIONAL, true_term);

  // Commit accumulated wire definitions via macro substitution.
  // Aliased entries belong in the parent's scope; everything else
  // uses the current prefix.
  for (auto & [sym, term] : pending_comb_updates_) {
    string name;
    if (pending_comb_aliased_.count(sym)) {
      name = parent_prefix_.empty() ? string(sym->name)
                                    : parent_prefix_ + "." + string(sym->name);
    } else {
      name = make_name(string(sym->name));
    }
    symbol_to_term_[sym] = term;
    fts_.name_term(name, term);
    logger.log(2, "SystemVerilogEncoder: always_comb (wire) {} := ...", name);
  }
  pending_comb_updates_.clear();
  pending_comb_aliased_.clear();
}

void SystemVerilogEncoder::process_initial(
    const slang::ast::ProceduralBlockSymbol & proc)
{
  Term true_term = solver_->make_term(true);
  process_statement(proc.getBody(), StmtContext::INITIAL, true_term);
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

  auto desc = resolve_lvalue(lhs_expr, eval_ctx());
  if (!desc) return;
  const Symbol * sym = desc->base;
  auto alias_it = port_output_aliases_.find(sym);
  bool aliased = alias_it != port_output_aliases_.end();
  if (aliased) sym = alias_it->second;

  uint64_t lo = desc->lo;
  uint64_t hi = desc->hi;
  uint64_t slice_w = hi - lo + 1;

  Term rhs = expr_to_term(rhs_expr);
  rhs = resize_to(rhs, slice_w);

  // Wire LHS: macro-substitute the *full-width* defining expression.
  // For a partial LHS (`assign arr[i] = ...`) we splice the slice
  // into whatever was previously stored under `sym`; for the very
  // first partial assign there is no prior term to splice into, so
  // we fall through to the add_invar path below where the base has
  // already been declared as a free input.
  if (wire_symbols_.count(sym)) {
    bool full_write = (lo == 0 && hi + 1 == slice_w);
    Term new_term;
    if (full_write) {
      new_term = rhs;
    } else {
      auto sit = symbol_to_term_.find(sym);
      if (sit == symbol_to_term_.end()) {
        // No seed term: partial assigns to a wire base aren't
        // expressible without it, so this base shouldn't have been
        // classified as a wire.  Skip.
        return;
      }
      new_term = replace_bits(sit->second, rhs, lo, hi);
    }
    string name;
    if (aliased) {
      name = parent_prefix_.empty() ? string(sym->name)
                                    : parent_prefix_ + "." + string(sym->name);
    } else {
      name = make_name(string(sym->name));
    }
    symbol_to_term_[sym] = new_term;
    fts_.name_term(name, new_term);
    logger.log(
        2, "SystemVerilogEncoder: continuous assign (wire) {} := ...", name);
    return;
  }

  // Fallback: existing variable (e.g., output port reg, or a
  // partially-driven base that wasn't classified as a wire).
  // Constrain the appropriate slice via add_constraint (which
  // tolerates input vars in the term).
  auto it = symbol_to_term_.find(sym);
  if (it != symbol_to_term_.end()) {
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

void SystemVerilogEncoder::process_instance(
    const slang::ast::InstanceSymbol & inst)
{
  using namespace slang::ast;

  // Switch the naming context so any state vars / wires declared
  // inside the child get hierarchical names like "<parent>.<inst>.<x>".
  // Also track the *parent's* prefix so wires redirected via
  // port_output_aliases_ (which live in the parent's scope) get
  // named correctly.
  string saved_prefix = prefix_;
  string saved_parent_prefix = parent_prefix_;
  parent_prefix_ = prefix_;
  prefix_ = saved_prefix + "." + string(inst.name);

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
      // left-hand side is the parent-side expression.
      if (conn_expr->kind == ExpressionKind::Assignment) {
        conn_expr = &conn_expr->as<AssignmentExpression>().left();
      }
      if (conn_expr->kind == ExpressionKind::NamedValue) {
        auto * parent_sym = &conn_expr->as<NamedValueExpression>().symbol;
        port_output_aliases_[internal] = parent_sym;
        output_aliases_added.push_back(internal);
      }
    } else {
      Term term = expr_to_term(*conn_expr);
      term = resize_to(term, port.getType().getBitWidth());
      symbol_to_term_[internal] = term;
      input_terms_added.push_back(internal);
    }
  }

  // Pre-scan the child's procedural blocks so internal NB-assigned
  // registers and blocking-assigned wires get classified before
  // declare_variables_internal runs.
  walk_members(inst.body, [&](const Symbol & m) {
    if (m.kind != SymbolKind::ProceduralBlock) return;
    auto & proc = m.as<ProceduralBlockSymbol>();
    if (proc.procedureKind == ProceduralBlockKind::AlwaysFF) {
      pre_scan_always_ff(proc.getBody());
    } else if (proc.procedureKind == ProceduralBlockKind::AlwaysComb) {
      pre_scan_always_comb(proc.getBody());
    } else if (proc.procedureKind == ProceduralBlockKind::Always) {
      std::unordered_set<const Symbol *> nb_targets;
      collect_nonblocking_targets(proc.getBody(), nb_targets);
      if (!nb_targets.empty()) {
        pre_scan_always_ff(proc.getBody());
      } else {
        pre_scan_always_comb(proc.getBody());
      }
    }
  });

  // Declare the child's internal (non-port) variables with the new
  // hierarchical prefix; ports are already bound through the
  // connection map above.
  declare_variables_internal(inst.body);

  // Combinational pass over child's body (and any sub-instances).
  walk_members(inst.body, [&](const Symbol & m) {
    if (m.kind == SymbolKind::ContinuousAssign) {
      process_continuous_assign(m.as<ContinuousAssignSymbol>());
    } else if (m.kind == SymbolKind::ProceduralBlock) {
      auto & proc = m.as<ProceduralBlockSymbol>();
      if (proc.procedureKind == ProceduralBlockKind::AlwaysComb) {
        process_always_comb(proc);
      } else if (proc.procedureKind == ProceduralBlockKind::Always) {
        std::unordered_set<const Symbol *> targets;
        collect_nonblocking_targets(proc.getBody(), targets);
        if (targets.empty()) {
          process_always_comb(proc);
        }
      }
    } else if (m.kind == SymbolKind::Instance) {
      process_instance(m.as<InstanceSymbol>());
    }
  });

  // Sequential / initial pass.
  walk_members(inst.body, [&](const Symbol & m) {
    if (m.kind != SymbolKind::ProceduralBlock) return;
    auto & proc = m.as<ProceduralBlockSymbol>();
    switch (proc.procedureKind) {
      case ProceduralBlockKind::AlwaysFF: process_always_ff(proc); break;
      case ProceduralBlockKind::Initial: process_initial(proc); break;
      case ProceduralBlockKind::Always: {
        std::unordered_set<const Symbol *> targets;
        collect_nonblocking_targets(proc.getBody(), targets);
        if (!targets.empty()) process_always_ff(proc);
        break;
      }
      default: break;
    }
  });

  // Restore context: undo the per-instance bindings so that a sibling
  // (or repeated) instantiation of the same module can be processed
  // cleanly.
  for (auto * sym : output_aliases_added) port_output_aliases_.erase(sym);
  for (auto * sym : input_terms_added) symbol_to_term_.erase(sym);
  prefix_ = saved_prefix;
  parent_prefix_ = saved_parent_prefix;
}

// ============================================================================
// Statement processing
// ============================================================================

void SystemVerilogEncoder::process_statement(const slang::ast::Statement & stmt,
                                             StmtContext ctx,
                                             const Term & condition)
{
  using namespace slang::ast;

  switch (stmt.kind) {
    case StatementKind::ExpressionStatement: {
      auto & es = stmt.as<ExpressionStatement>();
      auto & expr = es.expr;

      if (expr.kind == ExpressionKind::Assignment) {
        auto & assign = expr.as<AssignmentExpression>();
        auto & lhs_expr = assign.left();
        auto & rhs_expr = assign.right();

        auto desc = resolve_lvalue(lhs_expr, eval_ctx());
        if (!desc) break;
        const Symbol * sym = desc->base;
        auto alias_it = port_output_aliases_.find(sym);
        bool aliased = alias_it != port_output_aliases_.end();
        if (aliased) sym = alias_it->second;

        uint64_t lo = desc->lo;
        uint64_t hi = desc->hi;
        uint64_t slice_w = hi - lo + 1;

        // Figure out the full-base "previous" term so that
        //  (a) compound assignments can read the slice via the
        //      LValueReference stash, and
        //  (b) partial writes can be composed via replace_bits.
        Term prev_base;
        bool wire_comb =
            ctx == StmtContext::COMBINATIONAL && wire_symbols_.count(sym);
        Term state_term;  // only valid when ctx == NEXT_STATE
        if (wire_comb) {
          auto pit = pending_comb_updates_.find(sym);
          if (pit != pending_comb_updates_.end()) prev_base = pit->second;
        } else if (ctx == StmtContext::NEXT_STATE) {
          auto sit = symbol_to_term_.find(sym);
          if (sit == symbol_to_term_.end()) break;
          state_term = sit->second;
          auto pit = pending_next_updates_.find(state_term);
          prev_base =
              (pit != pending_next_updates_.end()) ? pit->second : state_term;
        } else {
          // COMBINATIONAL non-wire or INITIAL: prev_base is the
          // current (constant) value of the LHS used only for
          // compound-assignment self-reference.
          auto sit = symbol_to_term_.find(sym);
          if (sit != symbol_to_term_.end()) prev_base = sit->second;
        }

        // Stash the slice value for any LValueReference inside rhs.
        Term saved_lvalue = current_lvalue_term_;
        if (prev_base) {
          uint64_t pw = prev_base->get_sort()->get_width();
          if (lo == 0 && hi == pw - 1) {
            current_lvalue_term_ = prev_base;
          } else {
            current_lvalue_term_ =
                solver_->make_term(Op(Extract, hi, lo), prev_base);
          }
        } else {
          current_lvalue_term_ = Term();
        }

        Term rhs = expr_to_term(rhs_expr);
        current_lvalue_term_ = saved_lvalue;

        rhs = resize_to(rhs, slice_w);

        if (wire_comb) {
          if (aliased) pending_comb_aliased_.insert(sym);
          // Compose new full-base value from prev_base + slice rhs.
          // If there's no prev_base yet, treat the slice as the new
          // value (matches the existing latch-free default for
          // first-time conditional writes).
          Term combined;
          if (prev_base) {
            combined = replace_bits(prev_base, rhs, lo, hi);
          } else if (lo == 0 && hi + 1 == slice_w) {
            combined = rhs;
          } else {
            // Partial write to a wire with no seed term: skip; the
            // unwritten bits would be undefined.
            logger.log(1,
                       "SystemVerilogEncoder: partial write to undriven "
                       "wire {} -- skipped",
                       string(sym->name));
            break;
          }
          if (condition == solver_->make_term(true)) {
            pending_comb_updates_[sym] = combined;
          } else {
            Term def = prev_base ? prev_base : combined;
            pending_comb_updates_[sym] =
                solver_->make_term(Ite, condition, combined, def);
          }
          break;
        }

        auto it = symbol_to_term_.find(sym);
        if (it == symbol_to_term_.end()) break;
        Term lhs_term = it->second;
        uint64_t base_w = lhs_term->get_sort()->get_width();
        bool full_write = (lo == 0 && hi == base_w - 1);

        switch (ctx) {
          case StmtContext::NEXT_STATE: {
            Term combined =
                full_write ? rhs : replace_bits(prev_base, rhs, lo, hi);
            Term update;
            if (condition == solver_->make_term(true)) {
              update = combined;
            } else {
              update = solver_->make_term(Ite, condition, combined, prev_base);
            }
            pending_next_updates_[state_term] = update;
            break;
          }
          case StmtContext::COMBINATIONAL: {
            // Non-wire LHS (e.g. output port reg, or a partially-
            // written base): constrain the appropriate slice via
            // add_constraint, which accepts terms involving input
            // vars (so RHSes that reference input ports work).
            Term lhs_slice =
                full_write ? lhs_term
                           : solver_->make_term(Op(Extract, hi, lo), lhs_term);
            Term eq = solver_->make_term(Equal, lhs_slice, rhs);
            if (condition != solver_->make_term(true)) {
              eq = solver_->make_term(Implies, condition, eq);
            }
            fts_.add_constraint(eq);
            break;
          }
          case StmtContext::INITIAL: {
            Term lhs_slice =
                full_write ? lhs_term
                           : solver_->make_term(Op(Extract, hi, lo), lhs_term);
            Term eq = solver_->make_term(Equal, lhs_slice, rhs);
            fts_.constrain_init(eq);
            break;
          }
        }
      }
      break;
    }

    case StatementKind::Block: {
      auto & block = stmt.as<BlockStatement>();
      for_each_stmt_in_block(block, [&](const Statement & s) {
        process_statement(s, ctx, condition);
      });
      break;
    }

    case StatementKind::Conditional: {
      auto & cond_stmt = stmt.as<ConditionalStatement>();

      // Get the condition expression.
      // ConditionalStatement has conditions span; for simple if, there
      // is one condition.
      Term cond_term = expr_to_term(*cond_stmt.conditions[0].expr);

      // Ensure condition is a single-bit BV (for use with Ite).
      uint64_t cond_width = cond_term->get_sort()->get_width();
      if (cond_width > 1) {
        // Non-zero reduction: cond != 0
        Term zero = solver_->make_term(0, cond_term->get_sort());
        cond_term = solver_->make_term(Distinct, cond_term, zero);
        cond_term = solver_->make_term(
            Ite,
            cond_term,
            solver_->make_term(1, solver_->make_sort(BV, 1)),
            solver_->make_term(0, solver_->make_sort(BV, 1)));
      }

      // Build then-condition and else-condition.
      Sort bv1 = solver_->make_sort(BV, 1);
      Term one = solver_->make_term(1, bv1);
      Term zero = solver_->make_term(0, bv1);

      Term then_cond;
      Term else_cond;
      if (condition == solver_->make_term(true)) {
        // If the outer condition is trivially true, the condition is
        // just the if-expression.
        then_cond = solver_->make_term(Equal, cond_term, one);
        else_cond = solver_->make_term(Equal, cond_term, zero);
      } else {
        Term cond_eq_one = solver_->make_term(Equal, cond_term, one);
        Term cond_eq_zero = solver_->make_term(Equal, cond_term, zero);
        then_cond = solver_->make_term(And, condition, cond_eq_one);
        else_cond = solver_->make_term(And, condition, cond_eq_zero);
      }

      process_statement(cond_stmt.ifTrue, ctx, then_cond);
      if (cond_stmt.ifFalse) {
        process_statement(*cond_stmt.ifFalse, ctx, else_cond);
      }
      break;
    }

    case StatementKind::Case: {
      auto & case_stmt = stmt.as<CaseStatement>();
      Term sel = expr_to_term(case_stmt.expr);

      for (auto & item : case_stmt.items) {
        // Build OR of all patterns matching this item.
        Term item_cond;
        for (auto expr : item.expressions) {
          Term pat = expr_to_term(*expr);
          pat = resize_to(pat, sel->get_sort()->get_width());
          Term match = solver_->make_term(Equal, sel, pat);
          item_cond =
              item_cond ? solver_->make_term(Or, item_cond, match) : match;
        }
        Term full_cond = (condition == solver_->make_term(true))
                             ? item_cond
                             : solver_->make_term(And, condition, item_cond);
        process_statement(*item.stmt, ctx, full_cond);
      }
      if (case_stmt.defaultCase) {
        // Default: none of the other cases match.
        // For simplicity, we just use the outer condition.
        // A complete implementation would negate all other item conditions.
        process_statement(*case_stmt.defaultCase, ctx, condition);
      }
      break;
    }

    case StatementKind::Timed: {
      // Skip timing control (e.g., @(posedge clk)) and process the body.
      auto & timed = stmt.as<TimedStatement>();
      process_statement(timed.stmt, ctx, condition);
      break;
    }

    case StatementKind::ConcurrentAssertion: {
      auto & ca = stmt.as<ConcurrentAssertionStatement>();
      // Only handle 'assert' (not 'assume', 'cover', etc.).
      if (ca.assertionKind == AssertionKind::Assert) {
        // Extract the property expression.  Unwrap clocking wrappers
        // to reach the underlying SimpleAssertionExpr.
        const AssertionExpr * ae = &ca.propertySpec;
        // Strip ClockingAssertionExpr wrapper if present.
        if (ae->kind == AssertionExprKind::Clocking) {
          ae = &ae->as<ClockingAssertionExpr>().expr;
        }
        if (ae->kind == AssertionExprKind::Simple) {
          auto & simple = ae->as<SimpleAssertionExpr>();
          Term prop = expr_to_term(simple.expr);
          // Convert BV1 result to a Boolean equality check.
          uint64_t pw = prop->get_sort()->get_width();
          if (pw == 1) {
            prop = solver_->make_term(
                Equal, prop, solver_->make_term(1, solver_->make_sort(BV, 1)));
          } else {
            // Non-zero means true.
            prop = solver_->make_term(
                Distinct, prop, solver_->make_term(0, prop->get_sort()));
          }
          propvec_.push_back(prop);
          logger.log(
              1,
              "SystemVerilogEncoder: extracted concurrent assertion property");
        } else {
          logger.log(1,
                     "SystemVerilogEncoder: skipping unsupported assertion "
                     "expression kind {}",
                     static_cast<int>(ae->kind));
        }
      }
      break;
    }

    case StatementKind::VariableDeclaration: {
      // Procedural local variable (`int x = ...`).  Treated as an
      // immutable per-block constant: evaluate the initializer once
      // and bind it in the slang EvalContext and our SMT-side
      // loop_var_terms_ map.  Any later procedural write (outside a
      // for-loop step) is out of scope and effectively ignored.
      auto & vds = stmt.as<VariableDeclStatement>();
      auto & sym = vds.symbol;
      slang::ConstantValue cv;
      if (auto * init = sym.getInitializer()) {
        cv = init->eval(eval_ctx());
      }
      if (cv.bad()) {
        cv = sym.getType().getDefaultValue();
      }
      if (!cv.isInteger()) {
        throw PonoException("SystemVerilogEncoder: non-integer local '"
                            + string(sym.name) + "'");
      }
      eval_ctx().createLocal(&sym, cv);
      auto svint = cv.integer();
      uint64_t width = sym.getType().getBitWidth();
      if (width == 0) width = svint.getBitWidth();
      if (width == 0) width = 32;
      Sort sort = solver_->make_sort(BV, width);
      svint.setSigned(false);
      string val_str = svint.toString(slang::LiteralBase::Decimal, false);
      loop_var_terms_[&sym] = solver_->make_term(val_str, sort, 10);
      break;
    }

    case StatementKind::ForLoop: {
      // Compile-time unroll the loop.  Slang has already validated
      // that the bounds and steps are constant expressions for the
      // synthesizable subset.
      auto & loop = stmt.as<ForLoopStatement>();
      std::vector<const ValueSymbol *> declared;

      auto bind_var = [&](const VariableSymbol & lv) {
        slang::ConstantValue cv;
        if (auto * init = lv.getInitializer()) {
          cv = init->eval(eval_ctx());
        }
        if (cv.bad()) {
          cv = lv.getType().getDefaultValue();
        }
        if (!cv.isInteger()) {
          throw PonoException("SystemVerilogEncoder: non-integer for-loop var '"
                              + string(lv.name) + "'");
        }
        eval_ctx().createLocal(&lv, cv);
        declared.push_back(&lv);
      };

      auto refresh_bv = [&](const VariableSymbol & lv) {
        auto * cur = eval_ctx().findLocal(&lv);
        if (!cur || !cur->isInteger()) {
          throw PonoException("SystemVerilogEncoder: for-loop var '"
                              + string(lv.name) + "' lost its constant value");
        }
        auto svint = cur->integer();
        uint64_t width = lv.getType().getBitWidth();
        if (width == 0) width = svint.getBitWidth();
        if (width == 0) width = 32;
        Sort sort = solver_->make_sort(BV, width);
        svint.setSigned(false);
        string val_str = svint.toString(slang::LiteralBase::Decimal, false);
        loop_var_terms_[&lv] = solver_->make_term(val_str, sort, 10);
      };

      for (auto * lv : loop.loopVars) bind_var(*lv);
      for (auto * init : loop.initializers) {
        if (init->eval(eval_ctx()).bad()) {
          throw PonoException(
              "SystemVerilogEncoder: for-loop initializer eval failed");
        }
      }

      constexpr size_t MAX_ITERS = 65536;
      for (size_t it = 0;; ++it) {
        if (it >= MAX_ITERS) {
          throw PonoException("SystemVerilogEncoder: for-loop exceeded "
                              + std::to_string(MAX_ITERS) + " iterations");
        }
        if (loop.stopExpr) {
          auto sv = loop.stopExpr->eval(eval_ctx());
          if (sv.bad()) {
            throw PonoException(
                "SystemVerilogEncoder: for-loop stop eval failed");
          }
          if (!sv.isTrue()) break;
        }
        for (auto * lv : loop.loopVars) refresh_bv(*lv);
        process_statement(loop.body, ctx, condition);
        for (auto * step : loop.steps) {
          if (step->eval(eval_ctx()).bad()) {
            throw PonoException(
                "SystemVerilogEncoder: for-loop step eval failed");
          }
        }
      }

      for (auto * sym : declared) {
        loop_var_terms_.erase(sym);
        eval_ctx().deleteLocal(sym);
      }
      break;
    }

    default:
      // Other statement kinds (loops, etc.): not supported in
      // synthesizable subset. Log a warning and skip.
      logger.log(1,
                 "SystemVerilogEncoder: skipping unsupported statement kind {}",
                 static_cast<int>(stmt.kind));
      break;
  }
}

// ============================================================================
// Expression conversion
// ============================================================================

Term SystemVerilogEncoder::expr_to_term(const slang::ast::Expression & expr)
{
  using namespace slang::ast;

  switch (expr.kind) {
    case ExpressionKind::NamedValue: {
      auto & nv = expr.as<NamedValueExpression>();
      return lookup_symbol(&nv.symbol);
    }

    case ExpressionKind::HierarchicalValue: {
      // Cross-scope dotted read (e.g. `child_inst.reg`).  Slang has
      // already resolved the dotted path to the target ValueSymbol;
      // lookup_symbol finds its term in the appropriate scope
      // provided the referenced instance has been encoded already.
      auto & hv = expr.as<HierarchicalValueExpression>();
      return lookup_symbol(&hv.symbol);
    }

    case ExpressionKind::LValueReference: {
      // Implicit self-reference produced by compound assignments
      // (`x &= y` -> `x = LValueReference & y`).  The owning
      // assignment handler must have stashed the current LHS term
      // before recursing into the RHS.
      if (!current_lvalue_term_) {
        throw PonoException(
            "SystemVerilogEncoder: LValueReference outside compound "
            "assignment");
      }
      return current_lvalue_term_;
    }

    case ExpressionKind::IntegerLiteral: {
      auto & lit = expr.as<IntegerLiteral>();
      uint64_t width = expr.type->getBitWidth();
      if (width == 0) width = 32;  // Default integer width.
      Sort sort = solver_->make_sort(BV, width);
      // Reinterpret the value as unsigned so that toString emits the raw
      // two's-complement bit pattern as a positive decimal.  Without
      // setSigned(false), signed-negative values would stringify as
      // "-N", which smt-switch's base-10 parser rejects.
      auto val = lit.getValue();
      val.setSigned(false);
      string val_str =
          val.toString(slang::LiteralBase::Decimal, /*includeBase=*/false);
      return solver_->make_term(val_str, sort, 10);
    }

    case ExpressionKind::UnbasedUnsizedIntegerLiteral: {
      auto & lit = expr.as<UnbasedUnsizedIntegerLiteral>();
      uint64_t width = expr.type->getBitWidth();
      if (width == 0) width = 1;
      Sort sort = solver_->make_sort(BV, width);
      auto val = lit.getValue();
      val.setSigned(false);
      string val_str =
          val.toString(slang::LiteralBase::Decimal, /*includeBase=*/false);
      return solver_->make_term(val_str, sort, 10);
    }

    case ExpressionKind::BinaryOp: {
      auto & binop = expr.as<BinaryExpression>();
      Term left = expr_to_term(binop.left());
      Term right = expr_to_term(binop.right());

      // Ensure operands have the same width.
      uint64_t lw = left->get_sort()->get_width();
      uint64_t rw = right->get_sort()->get_width();
      uint64_t max_w = max(lw, rw);
      left = resize_to(left, max_w);
      right = resize_to(right, max_w);

      uint64_t result_width = expr.type->getBitWidth();
      Term result;

      switch (binop.op) {
        case BinaryOperator::Add:
          result = solver_->make_term(BVAdd, left, right);
          break;
        case BinaryOperator::Subtract:
          result = solver_->make_term(BVSub, left, right);
          break;
        case BinaryOperator::Multiply:
          result = solver_->make_term(BVMul, left, right);
          break;
        case BinaryOperator::Divide:
          result = solver_->make_term(BVUdiv, left, right);
          break;
        case BinaryOperator::Mod:
          result = solver_->make_term(BVUrem, left, right);
          break;
        case BinaryOperator::BinaryAnd:
          result = solver_->make_term(BVAnd, left, right);
          break;
        case BinaryOperator::BinaryOr:
          result = solver_->make_term(BVOr, left, right);
          break;
        case BinaryOperator::BinaryXor:
          result = solver_->make_term(BVXor, left, right);
          break;
        case BinaryOperator::BinaryXnor: {
          Term xor_t = solver_->make_term(BVXor, left, right);
          result = solver_->make_term(BVNot, xor_t);
          break;
        }
        case BinaryOperator::Equality: {
          Term eq = solver_->make_term(Equal, left, right);
          Sort bv1 = solver_->make_sort(BV, 1);
          result = solver_->make_term(
              Ite, eq, solver_->make_term(1, bv1), solver_->make_term(0, bv1));
          break;
        }
        case BinaryOperator::Inequality: {
          Term eq = solver_->make_term(Equal, left, right);
          Sort bv1 = solver_->make_sort(BV, 1);
          result = solver_->make_term(
              Ite, eq, solver_->make_term(0, bv1), solver_->make_term(1, bv1));
          break;
        }
        case BinaryOperator::LessThan: {
          Term lt = solver_->make_term(BVUlt, left, right);
          Sort bv1 = solver_->make_sort(BV, 1);
          result = solver_->make_term(
              Ite, lt, solver_->make_term(1, bv1), solver_->make_term(0, bv1));
          break;
        }
        case BinaryOperator::LessThanEqual: {
          Term le = solver_->make_term(BVUle, left, right);
          Sort bv1 = solver_->make_sort(BV, 1);
          result = solver_->make_term(
              Ite, le, solver_->make_term(1, bv1), solver_->make_term(0, bv1));
          break;
        }
        case BinaryOperator::GreaterThan: {
          Term gt = solver_->make_term(BVUlt, right, left);
          Sort bv1 = solver_->make_sort(BV, 1);
          result = solver_->make_term(
              Ite, gt, solver_->make_term(1, bv1), solver_->make_term(0, bv1));
          break;
        }
        case BinaryOperator::GreaterThanEqual: {
          Term ge = solver_->make_term(BVUle, right, left);
          Sort bv1 = solver_->make_sort(BV, 1);
          result = solver_->make_term(
              Ite, ge, solver_->make_term(1, bv1), solver_->make_term(0, bv1));
          break;
        }
        case BinaryOperator::LogicalAnd: {
          // Logical AND: both operands nonzero.
          Sort bv1 = solver_->make_sort(BV, 1);
          Term l_nz = solver_->make_term(
              Distinct, left, solver_->make_term(0, left->get_sort()));
          Term r_nz = solver_->make_term(
              Distinct, right, solver_->make_term(0, right->get_sort()));
          Term both = solver_->make_term(And, l_nz, r_nz);
          result = solver_->make_term(Ite,
                                      both,
                                      solver_->make_term(1, bv1),
                                      solver_->make_term(0, bv1));
          break;
        }
        case BinaryOperator::LogicalOr: {
          Sort bv1 = solver_->make_sort(BV, 1);
          Term l_nz = solver_->make_term(
              Distinct, left, solver_->make_term(0, left->get_sort()));
          Term r_nz = solver_->make_term(
              Distinct, right, solver_->make_term(0, right->get_sort()));
          Term either = solver_->make_term(Or, l_nz, r_nz);
          result = solver_->make_term(Ite,
                                      either,
                                      solver_->make_term(1, bv1),
                                      solver_->make_term(0, bv1));
          break;
        }
        case BinaryOperator::LogicalShiftLeft:
          result = solver_->make_term(BVShl, left, right);
          break;
        case BinaryOperator::LogicalShiftRight:
          result = solver_->make_term(BVLshr, left, right);
          break;
        case BinaryOperator::ArithmeticShiftLeft:
          result = solver_->make_term(BVShl, left, right);
          break;
        case BinaryOperator::ArithmeticShiftRight:
          result = solver_->make_term(BVAshr, left, right);
          break;
        default:
          throw PonoException(
              "SystemVerilogEncoder: unsupported binary operator "
              + to_string(static_cast<int>(binop.op)));
      }
      return resize_to(result, result_width);
    }

    case ExpressionKind::UnaryOp: {
      auto & unop = expr.as<UnaryExpression>();
      Term operand = expr_to_term(unop.operand());
      uint64_t result_width = expr.type->getBitWidth();

      Term result;
      switch (unop.op) {
        case UnaryOperator::BitwiseNot:
          result = solver_->make_term(BVNot, operand);
          break;
        case UnaryOperator::LogicalNot: {
          Sort bv1 = solver_->make_sort(BV, 1);
          Term is_zero = solver_->make_term(
              Equal, operand, solver_->make_term(0, operand->get_sort()));
          result = solver_->make_term(Ite,
                                      is_zero,
                                      solver_->make_term(1, bv1),
                                      solver_->make_term(0, bv1));
          break;
        }
        case UnaryOperator::Minus:
          result = solver_->make_term(BVNeg, operand);
          break;
        case UnaryOperator::BitwiseAnd: {
          // Reduction AND: result is 1 if all bits are 1.
          Sort bv1 = solver_->make_sort(BV, 1);
          Term all_ones = solver_->make_term(
              Equal,
              operand,
              solver_->make_term(string(operand->get_sort()->get_width(), '1'),
                                 operand->get_sort(),
                                 2));
          result = solver_->make_term(Ite,
                                      all_ones,
                                      solver_->make_term(1, bv1),
                                      solver_->make_term(0, bv1));
          break;
        }
        case UnaryOperator::BitwiseOr: {
          // Reduction OR: result is 1 if any bit is 1.
          Sort bv1 = solver_->make_sort(BV, 1);
          Term any_one = solver_->make_term(
              Distinct, operand, solver_->make_term(0, operand->get_sort()));
          result = solver_->make_term(Ite,
                                      any_one,
                                      solver_->make_term(1, bv1),
                                      solver_->make_term(0, bv1));
          break;
        }
        case UnaryOperator::BitwiseXor: {
          // Reduction XOR: parity of bits. For a BV of width n,
          // XOR all bits together.
          uint64_t w = operand->get_sort()->get_width();
          Sort bv1 = solver_->make_sort(BV, 1);
          result = solver_->make_term(Op(Extract, 0, 0), operand);
          for (uint64_t i = 1; i < w; i++) {
            Term bit = solver_->make_term(Op(Extract, i, i), operand);
            result = solver_->make_term(BVXor, result, bit);
          }
          break;
        }
        default:
          throw PonoException(
              "SystemVerilogEncoder: unsupported unary operator "
              + to_string(static_cast<int>(unop.op)));
      }
      return resize_to(result, result_width);
    }

    case ExpressionKind::Conversion: {
      auto & conv = expr.as<ConversionExpression>();
      Term inner = expr_to_term(conv.operand());
      uint64_t target_width = expr.type->getBitWidth();
      return resize_to(inner, target_width);
    }

    case ExpressionKind::Concatenation: {
      auto & concat = expr.as<ConcatenationExpression>();
      auto operands = concat.operands();
      if (operands.empty()) {
        throw PonoException("SystemVerilogEncoder: empty concatenation");
      }
      // Concatenate from MSB (first) to LSB (last).
      Term result = expr_to_term(*operands[0]);
      for (size_t i = 1; i < operands.size(); i++) {
        Term next = expr_to_term(*operands[i]);
        result = solver_->make_term(Concat, result, next);
      }
      return result;
    }

    case ExpressionKind::Replication: {
      auto & repl = expr.as<ReplicationExpression>();
      Term inner = expr_to_term(repl.concat());
      // The count should be a compile-time constant.
      auto count_cv = repl.count().getConstant();
      if (!count_cv) {
        throw PonoException(
            "SystemVerilogEncoder: non-constant replication count");
      }
      auto count_opt = count_cv->integer().as<uint32_t>();
      if (!count_opt) {
        throw PonoException("SystemVerilogEncoder: invalid replication count");
      }
      uint32_t count = *count_opt;
      if (count == 0) {
        throw PonoException("SystemVerilogEncoder: zero replication count");
      }
      Term result = inner;
      for (uint32_t i = 1; i < count; i++) {
        result = solver_->make_term(Concat, result, inner);
      }
      return result;
    }

    case ExpressionKind::ElementSelect: {
      // Element select: `a[i]`.  For packed arrays the element width
      // can be more than one bit; the constant-index case uses
      // Extract with the full slice, while the dynamic case is
      // restricted to single-bit elements where shift+extract works.
      auto & sel = expr.as<ElementSelectExpression>();
      Term val = expr_to_term(sel.value());
      auto & sel_expr = sel.selector();
      uint64_t elem_w = expr.type->getBitWidth();
      if (elem_w == 0) elem_w = 1;

      // Try to evaluate the selector as a constant -- including the
      // case where the index is a loop counter bound via eval_ctx.
      std::optional<uint64_t> idx_const;
      if (sel_expr.getConstant()) {
        idx_const = sel_expr.getConstant()->integer().as<uint64_t>();
      } else {
        auto cv = sel_expr.eval(eval_ctx());
        if (cv.isInteger()) idx_const = cv.integer().as<uint64_t>();
      }

      if (idx_const) {
        uint64_t low = *idx_const * elem_w;
        uint64_t high = low + elem_w - 1;
        return solver_->make_term(Op(Extract, high, low), val);
      }

      // Dynamic select: shift right by (idx * elem_w) bits, then
      // extract the bottom `elem_w` bits.
      uint64_t val_w = val->get_sort()->get_width();
      Sort val_sort = solver_->make_sort(BV, val_w);
      Term idx = expr_to_term(sel_expr);
      idx = resize_to(idx, val_w);
      Term shift_amount = idx;
      if (elem_w != 1) {
        Term elem_w_term = solver_->make_term(elem_w, val_sort);
        shift_amount = solver_->make_term(BVMul, idx, elem_w_term);
      }
      Term shifted = solver_->make_term(BVLshr, val, shift_amount);
      return solver_->make_term(Op(Extract, elem_w - 1, 0), shifted);
    }

    case ExpressionKind::RangeSelect: {
      // Range select: a[hi:lo]
      auto & sel = expr.as<RangeSelectExpression>();
      Term val = expr_to_term(sel.value());
      auto & left_expr = sel.left();
      auto & right_expr = sel.right();

      // Both bounds should be compile-time constants for synthesizable code.
      if (!left_expr.getConstant() || !right_expr.getConstant()) {
        throw PonoException(
            "SystemVerilogEncoder: non-constant range select bounds");
      }
      auto hi_opt = left_expr.getConstant()->integer().as<uint64_t>();
      auto lo_opt = right_expr.getConstant()->integer().as<uint64_t>();
      if (!hi_opt || !lo_opt) {
        throw PonoException(
            "SystemVerilogEncoder: invalid range select bounds");
      }
      uint64_t hi = *hi_opt;
      uint64_t lo = *lo_opt;
      if (hi < lo) swap(hi, lo);
      return solver_->make_term(Op(Extract, hi, lo), val);
    }

    case ExpressionKind::ConditionalOp: {
      auto & ternary = expr.as<ConditionalExpression>();
      Term cond = expr_to_term(*ternary.conditions[0].expr);
      Term then_val = expr_to_term(ternary.left());
      Term else_val = expr_to_term(ternary.right());

      // Ensure then/else have the same width.
      uint64_t tw = then_val->get_sort()->get_width();
      uint64_t ew = else_val->get_sort()->get_width();
      uint64_t max_w = max(tw, ew);
      then_val = resize_to(then_val, max_w);
      else_val = resize_to(else_val, max_w);

      // Convert condition to Bool for Ite.
      uint64_t cw = cond->get_sort()->get_width();
      Term bool_cond;
      if (cw == 1) {
        bool_cond = solver_->make_term(
            Equal, cond, solver_->make_term(1, solver_->make_sort(BV, 1)));
      } else {
        bool_cond = solver_->make_term(
            Distinct, cond, solver_->make_term(0, cond->get_sort()));
      }
      return solver_->make_term(Ite, bool_cond, then_val, else_val);
    }

    default:
      throw PonoException("SystemVerilogEncoder: unsupported expression kind "
                          + to_string(static_cast<int>(expr.kind)));
  }
}

// ============================================================================
// Type conversion
// ============================================================================

Sort SystemVerilogEncoder::type_to_sort(const slang::ast::Type & type)
{
  if (type.isIntegral()) {
    uint64_t width = type.getBitWidth();
    if (width == 0) {
      throw PonoException("SystemVerilogEncoder: zero-width integral type");
    }
    return solver_->make_sort(BV, width);
  }

  throw PonoException("SystemVerilogEncoder: unsupported type kind");
}

// ============================================================================
// Helpers
// ============================================================================

Term SystemVerilogEncoder::lookup_symbol(const slang::ast::Symbol * sym) const
{
  using namespace slang::ast;

  // Procedural for-loop counter: bound to a per-iteration BV
  // constant for the duration of the unrolling.
  auto lvt = loop_var_terms_.find(sym);
  if (lvt != loop_var_terms_.end()) {
    return lvt->second;
  }

  // If `sym` is a child instance's output-port internal, redirect to
  // the parent-side wire so reads resolve to its term.
  auto alias_it = port_output_aliases_.find(sym);
  if (alias_it != port_output_aliases_.end()) {
    sym = alias_it->second;
  }

  // Wire being defined in the enclosing always_comb block: return
  // the partial accumulated term so that read-modify-write patterns
  // (e.g. `popcount = popcount + din[i];` inside an unrolled for
  // loop) see the previously-written value.
  auto pending_it = pending_comb_updates_.find(sym);
  if (pending_it != pending_comb_updates_.end()) {
    return pending_it->second;
  }

  auto it = symbol_to_term_.find(sym);
  if (it != symbol_to_term_.end()) {
    return it->second;
  }

  // Parameter / localparam: slang has already evaluated the value at
  // elaboration time.  Materialize a fresh BV constant from it so
  // references to `REQUESTERS`, `MAX_COUNT`, etc. fold to literals.
  if (sym->kind == SymbolKind::Parameter) {
    auto & param = sym->as<ParameterSymbol>();
    const auto & cv = param.getValue();
    if (!cv.isInteger()) {
      throw PonoException("SystemVerilogEncoder: non-integer parameter '"
                          + string(sym->name) + "'");
    }
    auto val = cv.integer();
    uint64_t width = param.getType().getBitWidth();
    if (width == 0) width = val.getBitWidth();
    if (width == 0) width = 32;
    Sort sort = solver_->make_sort(BV, width);
    val.setSigned(false);
    string val_str =
        val.toString(slang::LiteralBase::Decimal, /*includeBase=*/false);
    return solver_->make_term(val_str, sort, 10);
  }

  throw PonoException("SystemVerilogEncoder: unknown symbol '"
                      + string(sym->name) + "'");
}

string SystemVerilogEncoder::make_name(const string & name) const
{
  if (prefix_.empty()) return name;
  return prefix_ + "." + name;
}

Term SystemVerilogEncoder::resize_to(const Term & t, uint64_t target_width)
{
  uint64_t current_width = t->get_sort()->get_width();
  if (current_width == target_width) {
    return t;
  }
  if (current_width < target_width) {
    // Zero-extend.
    return solver_->make_term(Op(Zero_Extend, target_width - current_width), t);
  }
  // Truncate (extract lower bits).
  return solver_->make_term(Op(Extract, target_width - 1, 0), t);
}

Term SystemVerilogEncoder::replace_bits(const Term & base,
                                        const Term & slice,
                                        uint64_t lo,
                                        uint64_t hi)
{
  uint64_t base_w = base->get_sort()->get_width();
  if (lo == 0 && hi == base_w - 1) return slice;
  std::vector<Term> parts;
  if (hi + 1 < base_w) {
    parts.push_back(solver_->make_term(Op(Extract, base_w - 1, hi + 1), base));
  }
  parts.push_back(slice);
  if (lo > 0) {
    parts.push_back(solver_->make_term(Op(Extract, lo - 1, 0), base));
  }
  Term result = parts[0];
  for (size_t i = 1; i < parts.size(); ++i) {
    result = solver_->make_term(Concat, result, parts[i]);
  }
  return result;
}

}  // namespace pono
