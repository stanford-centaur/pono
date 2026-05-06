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

#include "frontends/sv_encoder.h"

#include <cassert>
#include <iostream>

#include "slang/ast/ASTVisitor.h"
#include "slang/ast/Compilation.h"
#include "slang/ast/Expression.h"
#include "slang/ast/Symbol.h"
#include "slang/ast/expressions/AssignmentExpressions.h"
#include "slang/ast/expressions/LiteralExpressions.h"
#include "slang/ast/expressions/MiscExpressions.h"
#include "slang/ast/expressions/OperatorExpressions.h"
#include "slang/ast/expressions/SelectExpressions.h"
#include "slang/ast/Statement.h"
#include "slang/ast/statements/ConditionalStatements.h"
#include "slang/ast/statements/MiscStatements.h"
#include "slang/ast/expressions/AssertionExpr.h"
#include "slang/ast/symbols/BlockSymbols.h"
#include "slang/ast/symbols/CompilationUnitSymbols.h"
#include "slang/ast/symbols/InstanceSymbols.h"
#include "slang/ast/symbols/MemberSymbols.h"
#include "slang/ast/symbols/PortSymbols.h"
#include "slang/ast/symbols/VariableSymbols.h"
#include "slang/ast/types/AllTypes.h"
#include "slang/ast/types/Type.h"
#include "slang/diagnostics/DiagnosticEngine.h"
#include "slang/syntax/SyntaxTree.h"

#include "utils/logger.h"

using namespace smt;
using namespace std;

namespace pono {

// ============================================================================
// Construction / destruction
// ============================================================================

SVEncoder::SVEncoder(string filename, FunctionalTransitionSystem & fts)
    : fts_(fts), solver_(fts.solver())
{
  encode(filename);
}

SVEncoder::~SVEncoder() = default;

// ============================================================================
// Top-level encoding pipeline
// ============================================================================

void SVEncoder::encode(const string & filename)
{
  // Parse the SystemVerilog source file.
  // SyntaxTree::fromFile returns an expected<shared_ptr<SyntaxTree>, ...>.
  auto tree_result = slang::syntax::SyntaxTree::fromFile(filename);
  if (!tree_result) {
    throw PonoException("SVEncoder: failed to parse file: " + filename);
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
      diag_messages =
          slang::DiagnosticEngine::reportAll(*sm, diagnostics);
    } else {
      diag_messages = "(unable to format diagnostics)";
    }
    throw PonoException("SVEncoder: errors in SystemVerilog file:\n"
                        + diag_messages);
  }

  // Walk the top-level instances.
  auto & root = compilation_->getRoot();
  auto top_instances = root.topInstances;
  if (top_instances.empty()) {
    throw PonoException(
        "SVEncoder: no top-level module instances found in " + filename);
  }

  // Process the first top-level module.
  // Multi-top designs could be supported by iterating, but for model
  // checking we typically have a single top module.
  process_module(*top_instances[0]);
}

// ============================================================================
// Module processing
// ============================================================================

void SVEncoder::process_module(const slang::ast::InstanceSymbol & inst)
{
  string inst_name(inst.name);
  logger.log(1, "SVEncoder: processing module {}", inst_name);
  prefix_ = inst_name;

  auto & body = inst.body;

  // First pass: identify state variable symbols by scanning always_ff blocks
  // for non-blocking assignment targets, before declaring variables.
  for (auto & member : body.members()) {
    if (member.kind == slang::ast::SymbolKind::ProceduralBlock) {
      auto & proc =
          member.as<slang::ast::ProceduralBlockSymbol>();
      if (proc.procedureKind
              == slang::ast::ProceduralBlockKind::AlwaysFF
          || proc.procedureKind
                 == slang::ast::ProceduralBlockKind::Always) {
        // Pre-scan: find all non-blocking assignment targets.
        pre_scan_always_ff(proc.getBody());
      }
    }
  }

  // Second pass: declare all variables (ports + internal declarations).
  declare_variables(body);

  // Third pass: process behavioral code and continuous assignments.
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
          // The LHS of a non-blocking assignment is a state variable.
          auto & lhs = assign.left();
          if (lhs.kind == ExpressionKind::NamedValue) {
            auto & nv = lhs.as<NamedValueExpression>();
            targets.insert(&nv.symbol);
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
    default:
      // Other statement types: nothing to extract.
      break;
  }
}

}  // anonymous namespace

// This is a private helper called from process_module; it is not in the
// header to avoid exposing slang::ast::Statement there.
void SVEncoder::pre_scan_always_ff(const slang::ast::Statement & body)
{
  collect_nonblocking_targets(body, state_var_symbols_);
}

// ============================================================================
// Variable declaration (first pass)
// ============================================================================

void SVEncoder::declare_variables(
    const slang::ast::InstanceBodySymbol & body)
{
  using namespace slang::ast;

  // Process ports first.
  for (auto port_sym : body.getPortList()) {
    if (port_sym->kind == SymbolKind::Port) {
      process_port(port_sym->as<PortSymbol>());
    }
  }

  // Process internal variable declarations (non-port variables).
  for (auto & member : body.members()) {
    if (member.kind == SymbolKind::Variable) {
      auto & var = member.as<VariableSymbol>();
      // Skip if already declared via port processing.
      if (symbol_to_term_.count(&var)) continue;

      string name = make_name(string(var.name));
      Sort sort = type_to_sort(var.getType());

      if (state_var_symbols_.count(&var)) {
        // This is a register: create a state variable.
        Term sv = fts_.make_statevar(name, sort);
        symbol_to_term_[&var] = sv;
        fts_.name_term(name, sv);
        logger.log(2, "SVEncoder: state var {} : bv{}", name,
                   sort->get_width());
      } else {
        // This is a wire/combinational signal: create an input variable
        // as a placeholder. It will be constrained when we process
        // always_comb or continuous assign.
        Term iv = fts_.make_inputvar(name, sort);
        symbol_to_term_[&var] = iv;
        fts_.name_term(name, iv);
        logger.log(2, "SVEncoder: wire {} : bv{}", name,
                   sort->get_width());
      }
    } else if (member.kind == SymbolKind::Net) {
      auto & net = member.as<NetSymbol>();
      if (symbol_to_term_.count(&net)) continue;

      string name = make_name(string(net.name));
      Sort sort = type_to_sort(net.getType());

      Term iv = fts_.make_inputvar(name, sort);
      symbol_to_term_[&net] = iv;
      fts_.name_term(name, iv);
      logger.log(2, "SVEncoder: net {} : bv{}", name,
                 sort->get_width());
    }
  }
}

void SVEncoder::process_port(const slang::ast::PortSymbol & port)
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
      logger.log(2, "SVEncoder: input port {} : bv{}", name,
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
    logger.log(2, "SVEncoder: input port {} : bv{}", name,
               sort->get_width());
  } else {
    // Output or inout: check if it's a register.
    if (state_var_symbols_.count(internal)) {
      Term sv = fts_.make_statevar(name, sort);
      symbol_to_term_[internal] = sv;
      symbol_to_term_[&port] = sv;
      fts_.name_term(name, sv);
      logger.log(2, "SVEncoder: output port (reg) {} : bv{}", name,
                 sort->get_width());
    } else {
      Term iv = fts_.make_inputvar(name, sort);
      symbol_to_term_[internal] = iv;
      symbol_to_term_[&port] = iv;
      fts_.name_term(name, iv);
      logger.log(2, "SVEncoder: output port (wire) {} : bv{}", name,
                 sort->get_width());
    }
  }
}

// ============================================================================
// Assignment processing (second pass)
// ============================================================================

void SVEncoder::process_assignments(
    const slang::ast::InstanceBodySymbol & body)
{
  using namespace slang::ast;

  Term true_term = solver_->make_term(true);

  for (auto & member : body.members()) {
    if (member.kind == SymbolKind::ProceduralBlock) {
      auto & proc = member.as<ProceduralBlockSymbol>();
      switch (proc.procedureKind) {
        case ProceduralBlockKind::AlwaysFF:
          process_always_ff(proc);
          break;
        case ProceduralBlockKind::AlwaysComb:
          process_always_comb(proc);
          break;
        case ProceduralBlockKind::Initial:
          process_initial(proc);
          break;
        case ProceduralBlockKind::Always: {
          // Legacy 'always' blocks: infer whether this is sequential
          // or combinational based on whether the block contains
          // non-blocking assignments to state variables.
          std::unordered_set<const Symbol *> targets;
          collect_nonblocking_targets(proc.getBody(), targets);
          if (!targets.empty()) {
            process_always_ff(proc);
          } else {
            process_always_comb(proc);
          }
          break;
        }
        default:
          // AlwaysLatch, Final, etc. -- skip for now.
          break;
      }
    } else if (member.kind == SymbolKind::ContinuousAssign) {
      process_continuous_assign(member.as<ContinuousAssignSymbol>());
    }
  }
}

void SVEncoder::process_always_ff(
    const slang::ast::ProceduralBlockSymbol & proc)
{
  pending_next_updates_.clear();

  // Use a null condition to represent "unconditional".
  Term true_term = solver_->make_term(true);
  process_statement(proc.getBody(), StmtContext::NEXT_STATE, true_term);

  // Commit all pending next-state updates.
  for (auto & [state_term, next_expr] : pending_next_updates_) {
    fts_.assign_next(state_term, next_expr);
    logger.log(2, "SVEncoder: assign_next {} := ...",
               fts_.get_name(state_term));
  }
}

void SVEncoder::process_always_comb(
    const slang::ast::ProceduralBlockSymbol & proc)
{
  Term true_term = solver_->make_term(true);
  process_statement(proc.getBody(), StmtContext::COMBINATIONAL, true_term);
}

void SVEncoder::process_initial(
    const slang::ast::ProceduralBlockSymbol & proc)
{
  Term true_term = solver_->make_term(true);
  process_statement(proc.getBody(), StmtContext::INITIAL, true_term);
}

void SVEncoder::process_continuous_assign(
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

  Term rhs = expr_to_term(rhs_expr);

  // For continuous assigns, the LHS wire is constrained to equal the RHS.
  if (lhs_expr.kind == ExpressionKind::NamedValue) {
    auto & nv = lhs_expr.as<NamedValueExpression>();
    auto it = symbol_to_term_.find(&nv.symbol);
    if (it != symbol_to_term_.end()) {
      Term lhs_term = it->second;
      rhs = resize_to(rhs, lhs_term->get_sort()->get_width());
      // Constrain the wire to equal the RHS expression.
      Term eq = solver_->make_term(Equal, lhs_term, rhs);
      fts_.add_invar(eq);
      logger.log(2, "SVEncoder: continuous assign {} = ...",
                 fts_.get_name(lhs_term));
    }
  }
}

// ============================================================================
// Statement processing
// ============================================================================

void SVEncoder::process_statement(const slang::ast::Statement & stmt,
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
        Term rhs = expr_to_term(rhs_expr);

        if (lhs_expr.kind == ExpressionKind::NamedValue) {
          auto & nv = lhs_expr.as<NamedValueExpression>();
          auto it = symbol_to_term_.find(&nv.symbol);
          if (it == symbol_to_term_.end()) break;

          Term lhs_term = it->second;
          uint64_t lhs_width = lhs_term->get_sort()->get_width();
          rhs = resize_to(rhs, lhs_width);

          switch (ctx) {
            case StmtContext::NEXT_STATE: {
              // Build conditional next-state: if (condition) rhs else
              // current
              Term update;
              auto pit = pending_next_updates_.find(lhs_term);
              Term prev = (pit != pending_next_updates_.end())
                              ? pit->second
                              : lhs_term;  // default: hold current value
              // Check if condition is just 'true'
              if (condition == solver_->make_term(true)) {
                update = rhs;
              } else {
                update = solver_->make_term(Ite, condition, rhs, prev);
              }
              pending_next_updates_[lhs_term] = update;
              break;
            }
            case StmtContext::COMBINATIONAL: {
              // Constrain wire = rhs (under condition).
              Term eq = solver_->make_term(Equal, lhs_term, rhs);
              if (condition != solver_->make_term(true)) {
                eq = solver_->make_term(Implies, condition, eq);
              }
              fts_.add_invar(eq);
              break;
            }
            case StmtContext::INITIAL: {
              // Constrain initial state: state_var = value.
              Term eq = solver_->make_term(Equal, lhs_term, rhs);
              fts_.constrain_init(eq);
              break;
            }
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
            Ite, cond_term,
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
          item_cond = item_cond ? solver_->make_term(Or, item_cond, match)
                                : match;
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
                Equal, prop,
                solver_->make_term(1, solver_->make_sort(BV, 1)));
          } else {
            // Non-zero means true.
            prop = solver_->make_term(
                Distinct, prop,
                solver_->make_term(0, prop->get_sort()));
          }
          propvec_.push_back(prop);
          logger.log(1, "SVEncoder: extracted concurrent assertion property");
        } else {
          logger.log(
              1,
              "SVEncoder: skipping unsupported assertion expression kind {}",
              static_cast<int>(ae->kind));
        }
      }
      break;
    }

    default:
      // Other statement kinds (loops, etc.): not supported in
      // synthesizable subset. Log a warning and skip.
      logger.log(1, "SVEncoder: skipping unsupported statement kind {}",
                 static_cast<int>(stmt.kind));
      break;
  }
}

// ============================================================================
// Expression conversion
// ============================================================================

Term SVEncoder::expr_to_term(const slang::ast::Expression & expr)
{
  using namespace slang::ast;

  switch (expr.kind) {
    case ExpressionKind::NamedValue: {
      auto & nv = expr.as<NamedValueExpression>();
      return lookup_symbol(&nv.symbol);
    }

    case ExpressionKind::IntegerLiteral: {
      auto & lit = expr.as<IntegerLiteral>();
      uint64_t width = expr.type->getBitWidth();
      if (width == 0) width = 32;  // Default integer width.
      Sort sort = solver_->make_sort(BV, width);
      // SVInt::toString(base, includeBase) returns std::string.
      auto val = lit.getValue();
      string val_str = val.toString(10, /* includeBase */ false);
      // If the value is negative (signed), handle via two's complement.
      // smt-switch make_term with base 10 expects unsigned decimal.
      if (!val_str.empty() && val_str[0] == '-') {
        // Negative value: compute two's complement as a positive number.
        // Use the slang SVInt directly: get the raw bits as hex.
        val_str = val.toString(16, false);
        return solver_->make_term(val_str, sort, 16);
      }
      return solver_->make_term(val_str, sort, 10);
    }

    case ExpressionKind::UnbasedUnsizedIntegerLiteral: {
      auto & lit = expr.as<UnbasedUnsizedIntegerLiteral>();
      uint64_t width = expr.type->getBitWidth();
      if (width == 0) width = 1;
      Sort sort = solver_->make_sort(BV, width);
      auto val = lit.getValue();
      string val_str = val.toString(10, false);
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
              Ite, eq, solver_->make_term(1, bv1),
              solver_->make_term(0, bv1));
          break;
        }
        case BinaryOperator::Inequality: {
          Term eq = solver_->make_term(Equal, left, right);
          Sort bv1 = solver_->make_sort(BV, 1);
          result = solver_->make_term(
              Ite, eq, solver_->make_term(0, bv1),
              solver_->make_term(1, bv1));
          break;
        }
        case BinaryOperator::LessThan: {
          Term lt = solver_->make_term(BVUlt, left, right);
          Sort bv1 = solver_->make_sort(BV, 1);
          result = solver_->make_term(
              Ite, lt, solver_->make_term(1, bv1),
              solver_->make_term(0, bv1));
          break;
        }
        case BinaryOperator::LessThanEqual: {
          Term le = solver_->make_term(BVUle, left, right);
          Sort bv1 = solver_->make_sort(BV, 1);
          result = solver_->make_term(
              Ite, le, solver_->make_term(1, bv1),
              solver_->make_term(0, bv1));
          break;
        }
        case BinaryOperator::GreaterThan: {
          Term gt = solver_->make_term(BVUlt, right, left);
          Sort bv1 = solver_->make_sort(BV, 1);
          result = solver_->make_term(
              Ite, gt, solver_->make_term(1, bv1),
              solver_->make_term(0, bv1));
          break;
        }
        case BinaryOperator::GreaterThanEqual: {
          Term ge = solver_->make_term(BVUle, right, left);
          Sort bv1 = solver_->make_sort(BV, 1);
          result = solver_->make_term(
              Ite, ge, solver_->make_term(1, bv1),
              solver_->make_term(0, bv1));
          break;
        }
        case BinaryOperator::LogicalAnd: {
          // Logical AND: both operands nonzero.
          Sort bv1 = solver_->make_sort(BV, 1);
          Term l_nz = solver_->make_term(
              Distinct, left,
              solver_->make_term(0, left->get_sort()));
          Term r_nz = solver_->make_term(
              Distinct, right,
              solver_->make_term(0, right->get_sort()));
          Term both = solver_->make_term(And, l_nz, r_nz);
          result = solver_->make_term(
              Ite, both, solver_->make_term(1, bv1),
              solver_->make_term(0, bv1));
          break;
        }
        case BinaryOperator::LogicalOr: {
          Sort bv1 = solver_->make_sort(BV, 1);
          Term l_nz = solver_->make_term(
              Distinct, left,
              solver_->make_term(0, left->get_sort()));
          Term r_nz = solver_->make_term(
              Distinct, right,
              solver_->make_term(0, right->get_sort()));
          Term either = solver_->make_term(Or, l_nz, r_nz);
          result = solver_->make_term(
              Ite, either, solver_->make_term(1, bv1),
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
              "SVEncoder: unsupported binary operator "
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
              Equal, operand,
              solver_->make_term(0, operand->get_sort()));
          result = solver_->make_term(
              Ite, is_zero, solver_->make_term(1, bv1),
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
              Equal, operand,
              solver_->make_term(
                  string(operand->get_sort()->get_width(), '1'),
                  operand->get_sort(), 2));
          result = solver_->make_term(
              Ite, all_ones, solver_->make_term(1, bv1),
              solver_->make_term(0, bv1));
          break;
        }
        case UnaryOperator::BitwiseOr: {
          // Reduction OR: result is 1 if any bit is 1.
          Sort bv1 = solver_->make_sort(BV, 1);
          Term any_one = solver_->make_term(
              Distinct, operand,
              solver_->make_term(0, operand->get_sort()));
          result = solver_->make_term(
              Ite, any_one, solver_->make_term(1, bv1),
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
              "SVEncoder: unsupported unary operator "
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
        throw PonoException("SVEncoder: empty concatenation");
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
        throw PonoException("SVEncoder: non-constant replication count");
      }
      auto count_opt = count_cv->integer().as<uint32_t>();
      if (!count_opt) {
        throw PonoException("SVEncoder: invalid replication count");
      }
      uint32_t count = *count_opt;
      if (count == 0) {
        throw PonoException("SVEncoder: zero replication count");
      }
      Term result = inner;
      for (uint32_t i = 1; i < count; i++) {
        result = solver_->make_term(Concat, result, inner);
      }
      return result;
    }

    case ExpressionKind::ElementSelect: {
      // Single bit select: a[i]
      auto & sel = expr.as<ElementSelectExpression>();
      Term val = expr_to_term(sel.value());
      auto & sel_expr = sel.selector();

      // If the selector is a compile-time constant, use Extract.
      if (sel_expr.getConstant()) {
        auto idx_opt = sel_expr.getConstant()->integer().as<uint64_t>();
        if (idx_opt) {
          uint64_t idx = *idx_opt;
          return solver_->make_term(Op(Extract, idx, idx), val);
        }
      }

      // Dynamic bit select: shift right by index, then extract bit 0.
      Term idx = expr_to_term(sel_expr);
      idx = resize_to(idx, val->get_sort()->get_width());
      Term shifted = solver_->make_term(BVLshr, val, idx);
      return solver_->make_term(Op(Extract, 0, 0), shifted);
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
            "SVEncoder: non-constant range select bounds");
      }
      auto hi_opt = left_expr.getConstant()->integer().as<uint64_t>();
      auto lo_opt = right_expr.getConstant()->integer().as<uint64_t>();
      if (!hi_opt || !lo_opt) {
        throw PonoException(
            "SVEncoder: invalid range select bounds");
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
            Equal, cond,
            solver_->make_term(1, solver_->make_sort(BV, 1)));
      } else {
        bool_cond = solver_->make_term(
            Distinct, cond,
            solver_->make_term(0, cond->get_sort()));
      }
      return solver_->make_term(Ite, bool_cond, then_val, else_val);
    }

    default:
      throw PonoException(
          "SVEncoder: unsupported expression kind "
          + to_string(static_cast<int>(expr.kind)));
  }
}

// ============================================================================
// Type conversion
// ============================================================================

Sort SVEncoder::type_to_sort(const slang::ast::Type & type)
{
  if (type.isIntegral()) {
    uint64_t width = type.getBitWidth();
    if (width == 0) {
      throw PonoException("SVEncoder: zero-width integral type");
    }
    return solver_->make_sort(BV, width);
  }

  throw PonoException("SVEncoder: unsupported type kind");
}

// ============================================================================
// Helpers
// ============================================================================

Term SVEncoder::lookup_symbol(const slang::ast::Symbol * sym) const
{
  auto it = symbol_to_term_.find(sym);
  if (it != symbol_to_term_.end()) {
    return it->second;
  }
  throw PonoException("SVEncoder: unknown symbol '"
                      + string(sym->name) + "'");
}

string SVEncoder::make_name(const string & name) const
{
  if (prefix_.empty()) return name;
  return prefix_ + "." + name;
}

Term SVEncoder::resize_to(const Term & t, uint64_t target_width)
{
  uint64_t current_width = t->get_sort()->get_width();
  if (current_width == target_width) {
    return t;
  }
  if (current_width < target_width) {
    // Zero-extend.
    return solver_->make_term(
        Op(Zero_Extend, target_width - current_width), t);
  }
  // Truncate (extract lower bits).
  return solver_->make_term(Op(Extract, target_width - 1, 0), t);
}

}  // namespace pono
