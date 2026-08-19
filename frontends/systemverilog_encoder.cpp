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
 ** SVA design decisions
 ** ---------------------
 ** A few SVA constructs have no single "obviously correct" encoding given
 ** this encoder's model (a single global clock, and a propvec()-of-safety-
 ** properties / ltl_justice()-of-liveness-obligations interface with no
 ** notion of coverage or multiple clock domains). The choices made here are
 ** deliberate and documented so future changes don't accidentally drift
 ** from them:
 **
 **   - `cover property (P)` / immediate `cover (P)`: modeled via
 **     reachability duality -- checked exactly like `assert property (!P)`
 **     (or `assert (!P)`), so a "violation" of that surrogate assertion is
 **     precisely "P was reached". This is the standard way to expose a
 **     coverage goal through a safety-property-only interface. Temporal/
 **     sequence-shaped cover goals are out of scope (negating a liveness
 **     obligation isn't a reachability check) and throw a clear error.
 **
 **   - Multiclock properties (a property whose sequence mentions more than
 **     one `@(posedge clkN)`): this encoder has no clock-domain-crossing
 **     model at all -- every design in this frontend's test suite already
 **     implicitly assumes one global clock advances the whole design by one
 **     cycle per sample. Rather than reject a multiclock property outright,
 **     every named clock is treated identically (each clock edge mentioned
 **     anywhere in a sequence advances the same global pono-cycle). This is
 **     the minimal behavior consistent with the encoder's existing single-
 **     clock assumption, not a real multi-domain/clock-ratio model.
 **/

#include "frontends/systemverilog_encoder.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
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
#include "slang/ast/expressions/CallExpression.h"
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
#include "slang/ast/types/AllTypes.h"
#include "slang/ast/types/Type.h"
#include "slang/diagnostics/DiagnosticEngine.h"
#include "slang/numeric/SVInt.h"
#include "slang/syntax/AllSyntax.h"
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

// Returns the source label of a concurrent assertion statement (e.g. the
// `p1` in `p1: assert property (...)`), or "<unnamed>" if it has none.
std::string assertion_label(const slang::ast::Statement & stmt);

// Internal control-flow signal for `break`/`continue`/`disable`,
// thrown by process_statement()'s Break/Continue/Disable cases and
// caught by whichever enclosing construct can absorb it: a ForLoop
// (for Break/Continue) or a named Block whose symbol matches
// disable_target (for Disable). This only correctly models
// compile-time-reachable control flow -- e.g. `if (i == 2) break;`
// where `i` is an already-unrolled `for`-loop counter, handled by the
// Conditional case's constant-fold fast path, which branches in C++
// rather than building a symbolic guard for both arms, so this signal
// can propagate out of only the taken branch exactly like a real
// break/continue/disable would.  A signal thrown from inside a
// *runtime*-dependent condition (which the general symbolic-guard
// path processes unconditionally for both arms) has no correct
// interpretation here and is never meant to be caught by anything;
// process_always_ff()/process_always_comb()/process_initial() convert
// any instance that escapes all the way out into a clear
// PonoException rather than letting it propagate as a raw internal
// exception type.
struct LoopControlSignal
{
  enum Kind
  {
    Break,
    Continue,
    Disable
  } kind;
  const slang::ast::Symbol * disable_target = nullptr;  // only for Disable
};
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
    } else if (m.kind == SymbolKind::InstanceArray) {
      // Arrayed instantiation (`mod inst[N-1:0] (...)`): slang
      // creates one child Instance per element, each with its own
      // port connections already sliced to the correct bus range
      // (see AssignmentExpressions.cpp's use of InstanceSymbol::
      // arrayPath). Flatten them into the same dispatch as a plain
      // Instance, pushing a bracket-indexed prefix so each element's
      // state gets a unique hierarchical name -- the elements
      // themselves are unnamed (only the array is named).
      auto & arr = m.as<InstanceArraySymbol>();
      std::string saved_prefix = prefix_;
      std::string arr_name = std::string(arr.name);
      for (size_t i = 0; i < arr.elements.size(); ++i) {
        auto * element = arr.elements[i];
        if (!element) continue;
        prefix_ = saved_prefix + "." + arr_name + "[" + std::to_string(i) + "]";
        fn(*element);
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

SystemVerilogEncoder::SystemVerilogEncoder(
    string filename,
    FunctionalTransitionSystem & fts,
    const std::vector<std::string> & filelists)
    : fts_(fts), solver_(fts.solver())
{
  encode(filename, filelists);
}

SystemVerilogEncoder::~SystemVerilogEncoder() = default;

// ============================================================================
// List-file (".f") parsing
// ============================================================================

namespace {

// Parses a SystemVerilog list ("dot-f") file: one source-file path per
// line, blank lines and '#'/"//" comment lines are skipped, and relative
// paths resolve against the directory containing `filelist_path`. Lines
// that look like a tool directive (starting with '+' or '-', e.g.
// "+incdir+" or "-y") are rejected rather than silently treated as
// filenames, since those directives are not supported.
std::vector<std::string> parse_dot_f_file(const std::string & filelist_path)
{
  std::ifstream in(filelist_path);
  if (!in) {
    throw PonoException("SystemVerilogEncoder: failed to open list file: "
                        + filelist_path);
  }

  std::filesystem::path base_dir =
      std::filesystem::path(filelist_path).parent_path();

  std::vector<std::string> files;
  std::string line;
  while (std::getline(in, line)) {
    size_t begin = line.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) {
      continue;  // blank line
    }
    size_t end = line.find_last_not_of(" \t\r\n");
    std::string trimmed = line.substr(begin, end - begin + 1);

    if (trimmed.rfind('#', 0) == 0 || trimmed.rfind("//", 0) == 0) {
      continue;  // comment line
    }
    if (trimmed[0] == '+' || trimmed[0] == '-') {
      throw PonoException(
          "SystemVerilogEncoder: unsupported list-file directive '" + trimmed
          + "' in " + filelist_path
          + " -- only bare source file paths are supported.");
    }

    std::filesystem::path p(trimmed);
    if (p.is_relative()) {
      p = base_dir / p;
    }
    files.push_back(p.string());
  }

  return files;
}

}  // namespace

// ============================================================================
// Top-level encoding pipeline
// ============================================================================

void SystemVerilogEncoder::encode(const string & filename,
                                  const std::vector<std::string> & filelists)
{
  // Build the full list of source files: the primary file followed by
  // every file named in each list-file, in order.
  std::vector<std::string> source_files{ filename };
  for (const auto & filelist : filelists) {
    std::vector<std::string> expanded = parse_dot_f_file(filelist);
    source_files.insert(source_files.end(), expanded.begin(), expanded.end());
  }

  // Create a compilation and add every source file's syntax tree to it, so
  // they are all elaborated together.
  compilation_ = make_unique<slang::ast::Compilation>();
  for (const auto & source_file : source_files) {
    // SyntaxTree::fromFile returns an expected<shared_ptr<SyntaxTree>, ...>.
    auto tree_result = slang::syntax::SyntaxTree::fromFile(source_file);
    if (!tree_result) {
      throw PonoException("SystemVerilogEncoder: failed to parse file: "
                          + source_file);
    }
    compilation_->addSyntaxTree(std::move(tree_result).value());
  }

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

  // First pass: identify state variable symbols by scanning always_ff
  // blocks for non-blocking assignment targets, before declaring any
  // variables anywhere in the design -- recurses into every descendant
  // instance up front (not just this body's own direct members) so a
  // sibling instance visited earlier in source order (e.g. an
  // `interface` instance whose members are actually driven by a later
  // sibling's always_ff through a hierarchical/interface-port
  // reference) doesn't get its members wrongly declared as free inputs
  // before its true driver is discovered.
  pre_scan_state_vars(body);

  // Second pre-pass: identify combinational wire symbols from
  // always_comb blocks, legacy `always` blocks without non-blocking
  // assignments, and continuous-assign LHS values.  Wires are not
  // independent variables -- they will be macro-substituted with their
  // defining expressions, so we must skip declaring them as input vars.
  walk_members(body, [&](const slang::ast::Symbol & member) {
    if (member.kind == slang::ast::SymbolKind::ProceduralBlock) {
      auto & proc = member.as<slang::ast::ProceduralBlockSymbol>();
      if (proc.procedureKind == slang::ast::ProceduralBlockKind::AlwaysComb) {
        pre_scan_always_comb(proc.getBody(), proc);
      } else if (proc.procedureKind
                 == slang::ast::ProceduralBlockKind::Always) {
        // Legacy always: combinational iff it has no non-blocking
        // assignments to identify it as sequential.
        std::unordered_set<const slang::ast::Symbol *> nb_targets;
        collect_nonblocking_targets(proc.getBody(), nb_targets);
        if (nb_targets.empty()) {
          pre_scan_always_comb(proc.getBody(), proc);
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
            wire_drivers_[sym] = { &ca, nullptr, prefix_, parent_prefix_ };
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
    case StatementKind::WhileLoop:
      collect_blocking_targets(
          stmt.as<WhileLoopStatement>().body, full, partial);
      break;
    case StatementKind::DoWhileLoop:
      collect_blocking_targets(
          stmt.as<DoWhileLoopStatement>().body, full, partial);
      break;
    case StatementKind::RepeatLoop:
      collect_blocking_targets(
          stmt.as<RepeatLoopStatement>().body, full, partial);
      break;
    case StatementKind::ForeachLoop:
      collect_blocking_targets(
          stmt.as<ForeachLoopStatement>().body, full, partial);
      break;
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
    case StatementKind::WhileLoop:
      collect_nonblocking_targets(stmt.as<WhileLoopStatement>().body, targets);
      break;
    case StatementKind::DoWhileLoop:
      collect_nonblocking_targets(stmt.as<DoWhileLoopStatement>().body,
                                  targets);
      break;
    case StatementKind::RepeatLoop:
      collect_nonblocking_targets(stmt.as<RepeatLoopStatement>().body, targets);
      break;
    case StatementKind::ForeachLoop:
      collect_nonblocking_targets(stmt.as<ForeachLoopStatement>().body,
                                  targets);
      break;
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
    case ExpressionKind::MemberAccess:
      return find_lhs_base(lhs.as<MemberAccessExpression>().value());
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
    case ExpressionKind::MemberAccess: {
      // Packed-struct/union field write (`s.field <= ...`): narrow the
      // inner base's range by the field's own bitOffset, mirroring the
      // read-side logic in expr_to_term()'s MemberAccess case. Additive
      // offsets compose correctly for nested access (`s.a.x`), since
      // each FieldSymbol's bitOffset is relative to its own immediately
      // enclosing struct/union type.
      auto & ma = lhs.as<MemberAccessExpression>();
      if (ma.member.kind != SymbolKind::Field) return std::nullopt;
      auto inner = resolve_lvalue(ma.value(), ctx);
      if (!inner) return std::nullopt;
      auto & field = ma.member.as<FieldSymbol>();
      uint64_t w = field.getType().getBitWidth();
      if (w == 0) return std::nullopt;
      uint64_t lo = inner->lo + field.bitOffset;
      uint64_t hi = lo + w - 1;
      if (hi > inner->hi) return std::nullopt;
      return LValueDesc{ inner->base, lo, hi, inner->base_w };
    }
    default: return std::nullopt;
  }
}

std::string assertion_label(const slang::ast::Statement & stmt)
{
  if (auto * syntax = stmt.syntax) {
    if (auto * ca_syntax =
            syntax
                ->as_if<slang::syntax::ConcurrentAssertionStatementSyntax>()) {
      if (ca_syntax->label) {
        return std::string(ca_syntax->label->name.valueText());
      }
    }
  }
  return "<unnamed>";
}

}  // anonymous namespace

// This is a private helper called from process_module; it is not in the
// header to avoid exposing slang::ast::Statement there.
void SystemVerilogEncoder::pre_scan_always_ff(
    const slang::ast::Statement & body)
{
  collect_nonblocking_targets(body, state_var_symbols_);
}

namespace {
// `initial forever @(...) body` is a legacy structural spelling of
// `always @(...) body` -- unlike a general `forever` (which has no
// static iteration bound at all and is a genuine architectural
// boundary of this encoder's compile-time-unrolling model), this
// specific shape runs its (timing-controlled) body exactly once per
// pono-cycle, exactly like an always_ff/always block does. Returns the
// forever loop's own (Timed) body if `stmt` matches this shape
// (allowing the single-statement Block wrapper slang gives an
// `initial` block's top-level statement), nullptr otherwise.
const slang::ast::Statement * as_forever_event_body(
    const slang::ast::Statement & stmt)
{
  using namespace slang::ast;
  const Statement * s = &stmt;
  if (s->kind == StatementKind::Block) {
    auto & block = s->as<BlockStatement>();
    auto & inner = block.body;
    if (inner.kind == StatementKind::List) {
      auto & list = inner.as<StatementList>();
      if (list.list.size() != 1) return nullptr;
      s = list.list[0];
    } else {
      s = &inner;
    }
  }
  if (s->kind != StatementKind::ForeverLoop) return nullptr;
  auto & forever_stmt = s->as<ForeverLoopStatement>();
  if (forever_stmt.body.kind != StatementKind::Timed) return nullptr;
  return &forever_stmt.body;
}
}  // namespace

void SystemVerilogEncoder::pre_scan_state_vars(
    const slang::ast::InstanceBodySymbol & body)
{
  using namespace slang::ast;
  walk_members(body, [&](const Symbol & member) {
    if (member.kind == SymbolKind::ProceduralBlock) {
      auto & proc = member.as<ProceduralBlockSymbol>();
      if (proc.procedureKind == ProceduralBlockKind::AlwaysFF
          || proc.procedureKind == ProceduralBlockKind::Always) {
        pre_scan_always_ff(proc.getBody());
      } else if (proc.procedureKind == ProceduralBlockKind::AlwaysLatch) {
        pre_scan_always_latch(proc.getBody());
      } else if (proc.procedureKind == ProceduralBlockKind::Initial) {
        if (auto * forever_body = as_forever_event_body(proc.getBody())) {
          pre_scan_always_ff(*forever_body);
        }
      }
    } else if (member.kind == SymbolKind::Instance) {
      pre_scan_state_vars(member.as<InstanceSymbol>().body);
    }
  });
}

void SystemVerilogEncoder::pre_scan_always_comb(
    const slang::ast::Statement & body,
    const slang::ast::ProceduralBlockSymbol & proc)
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
      wire_drivers_[sym] = { nullptr, &proc, prefix_, parent_prefix_ };
    }
  }
  for (auto * sym : partial) {
    if (!wire_symbols_.count(sym)) {
      state_var_symbols_.insert(sym);
    }
  }
}

void SystemVerilogEncoder::pre_scan_always_latch(
    const slang::ast::Statement & body)
{
  std::unordered_set<const slang::ast::Symbol *> full, partial;
  collect_blocking_targets(body, full, partial);
  for (auto * sym : full) state_var_symbols_.insert(sym);
  for (auto * sym : partial) state_var_symbols_.insert(sym);
}

void SystemVerilogEncoder::pre_scan_instance(
    const slang::ast::InstanceSymbol & inst)
{
  using namespace slang::ast;

  // Each output (or inout) port of the child is driven by the child's
  // logic.  In the parent's view the connected expression is usually a
  // simple NamedValue (e.g. `child c (.sum(y))`), but for an instance-
  // array element it's a constant-index slice of a parent-side bus
  // (e.g. `fifo_data_out[i]`, already resolved by slang) -- either way
  // the underlying base symbol becomes a wire.  More elaborate
  // connections (concatenations, non-constant indices, etc.) are not
  // yet supported.
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
    auto * parent_sym = find_lhs_base(*conn_expr);
    if (!parent_sym) continue;
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
      // skip declaring a separate term for it.  A register is the
      // exception: unlike a comb wire (whose term is filled in later
      // via macro substitution when its driving assignment is
      // processed), no later pass ever assigns a term to a bare
      // pass-through symbol, so a register's state var must be
      // created here, keyed under the fully resolved alias root.
      if (port_output_aliases_.count(&var)) {
        if (state_var_symbols_.count(&var)) {
          auto resolved = resolve_output_alias(&var);
          // A ranged alias (a register connected through a bus-
          // element instance-array connection) isn't supported here;
          // leave unconstrained rather than create a wrongly-sized
          // state var.
          if (!resolved.has_range && !symbol_to_term_.count(resolved.sym)) {
            const Symbol * root = resolved.sym;
            string name = make_name(string(var.name));
            Sort sort = type_to_sort(var.getType());
            Term sv = fts_.make_statevar(name, sort);
            symbol_to_term_[root] = sv;
            fts_.name_term(name, sv);
            logger.log(2,
                       "SystemVerilogEncoder: state var (aliased) {} : bv{}",
                       name,
                       sort->get_width());
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

// ============================================================================
// Assignment processing (second pass)
// ============================================================================

void SystemVerilogEncoder::process_assignments(
    const slang::ast::InstanceBodySymbol & body)
{
  using namespace slang::ast;

  // Track which module's scope we're in so assertions processed below
  // can look up that module's `default disable iff` (if any).  Saved
  // and restored around this call so that, after recursing into a
  // child instance below, subsequent siblings in the parent's own
  // member loop see the parent's scope again.
  const Scope * saved_scope = current_scope_;
  current_scope_ = &body;

  // Combinational definitions are processed first so that wires have a
  // term assigned in symbol_to_term_ before any always_ff or initial
  // block (or assertion) tries to reference them.  Child instances are
  // walked here too -- a child's continuous assigns / always_comb
  // blocks may drive parent-side wires that downstream parent code
  // references.
  walk_members(body, [&](const Symbol & member) {
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
    }
  });

  // Sequential and assertion-bearing blocks come second.
  walk_members(body, [&](const Symbol & member) {
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
  pending_next_updates_.clear();

  // Use a null condition to represent "unconditional".
  Term true_term = solver_->make_term(true);
  try {
    process_statement(body, StmtContext::NEXT_STATE, true_term);
  }
  catch (const LoopControlSignal &) {
    // Only a compile-time-constant break/continue/disable (caught by
    // an enclosing ForLoop/named Block) is supported; one reached via
    // a runtime-dependent condition, or with no matching enclosing
    // construct at all, escapes all the way here.
    throw PonoException(
        "SystemVerilogEncoder: break/continue/disable is only supported "
        "when its condition is a compile-time constant (e.g. depends only "
        "on already-unrolled for-loop counters)");
  }

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
  try {
    process_statement(proc.getBody(), StmtContext::COMBINATIONAL, true_term);
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
  try {
    process_statement(proc.getBody(), StmtContext::INITIAL, true_term);
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
  if (!processed_drivers_.insert(&proc).second) return;
  process_always_comb(proc);
}

void SystemVerilogEncoder::process_continuous_assign_once(
    const slang::ast::ContinuousAssignSymbol & ca)
{
  if (!processed_drivers_.insert(&ca).second) return;
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

  auto desc = resolve_lvalue(lhs_expr, eval_ctx());
  if (!desc) return;
  const Symbol * sym = desc->base;
  auto resolved = resolve_output_alias(sym);
  sym = resolved.sym;
  bool aliased = resolved.aliased;

  uint64_t lo = desc->lo;
  uint64_t hi = desc->hi;
  uint64_t slice_w = hi - lo + 1;
  if (resolved.has_range) {
    // `sym` (pre-alias) is exactly bits [resolved.lo, resolved.hi] of
    // the resolved root -- remap this write into the root's bits.
    lo += resolved.lo;
    hi += resolved.lo;
  }

  Term rhs = expr_to_term(rhs_expr);
  rhs = resize_to(rhs, slice_w);

  // Wire LHS: macro-substitute the *full-width* defining expression.
  // For a partial LHS (`assign arr[i] = ...`, or one element of an
  // instance array wired to a slice of a bus) we splice the slice
  // into whatever was previously stored under `sym`, creating a fresh
  // placeholder to splice into on the very first such write.
  if (wire_symbols_.count(sym)) {
    auto sit = symbol_to_term_.find(sym);
    bool full_write = !resolved.has_range && lo == 0
                      && (sit == symbol_to_term_.end()
                          || hi == sit->second->get_sort()->get_width() - 1);
    Term new_term;
    if (full_write) {
      new_term = rhs;
    } else {
      new_term = replace_bits(wire_seed_term(sym), rhs, lo, hi);
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
      auto desc = resolve_lvalue(*conn_expr, eval_ctx());
      if (desc) {
        bool full = (desc->lo == 0 && desc->hi + 1 == desc->base_w);
        port_output_aliases_[internal] =
            full ? OutputAliasTarget{ desc->base }
                 : OutputAliasTarget{ desc->base, true, desc->lo, desc->hi };
        output_aliases_added.push_back(internal);
      }
    } else {
      Term term = expr_to_term(*conn_expr);
      term = resize_to(term, port.getType().getBitWidth());
      symbol_to_term_[internal] = term;
      input_terms_added.push_back(internal);
    }
  }

  // Pre-scan the child's blocking-assigned wires so they get classified
  // before declare_variables_internal runs. NB-assigned registers are
  // *not* pre-scanned here -- process_module()'s pre_scan_state_vars()
  // already classified every always_ff/always block in the whole
  // design tree, including this instance's, before any instance's
  // variables were declared.
  walk_members(inst.body, [&](const Symbol & m) {
    if (m.kind != SymbolKind::ProceduralBlock) return;
    auto & proc = m.as<ProceduralBlockSymbol>();
    if (proc.procedureKind == ProceduralBlockKind::AlwaysComb) {
      pre_scan_always_comb(proc.getBody(), proc);
    } else if (proc.procedureKind == ProceduralBlockKind::Always) {
      std::unordered_set<const Symbol *> nb_targets;
      collect_nonblocking_targets(proc.getBody(), nb_targets);
      if (nb_targets.empty()) {
        pre_scan_always_comb(proc.getBody(), proc);
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
  walk_members(inst.body, [&](const Symbol & m) {
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
  for (auto * sym : output_aliases_added) port_output_aliases_.erase(sym);
  for (auto * sym : input_terms_added) symbol_to_term_.erase(sym);
  prefix_ = saved_prefix;
  parent_prefix_ = saved_parent_prefix;
}

// ============================================================================
// Statement processing
// ============================================================================

void SystemVerilogEncoder::process_dynamic_element_assign(
    const slang::ast::ElementSelectExpression & sel,
    const slang::ast::Expression & rhs_expr,
    StmtContext ctx,
    const Term & condition)
{
  using namespace slang::ast;

  // Only a direct `base[idx] = rhs` pattern is supported: the select
  // must sit directly on a plain variable, not on a nested select.
  if (sel.value().kind != ExpressionKind::NamedValue
      && sel.value().kind != ExpressionKind::HierarchicalValue) {
    return;
  }
  const Symbol * sym =
      (sel.value().kind == ExpressionKind::NamedValue)
          ? &sel.value().as<NamedValueExpression>().symbol
          : &sel.value().as<HierarchicalValueExpression>().symbol;
  auto resolved = resolve_output_alias(sym);
  sym = resolved.sym;
  bool aliased = resolved.aliased;
  if (resolved.has_range) {
    // Dynamic-index writes into a bus-element alias (e.g. one element
    // of an instance array wired to a slice of a parent bus) aren't
    // supported; leave unconstrained rather than risk a wrong
    // encoding.
    return;
  }

  uint64_t elem_w = sel.type->getBitWidth();
  if (elem_w == 0) elem_w = 1;

  bool wire_comb =
      ctx == StmtContext::COMBINATIONAL && wire_symbols_.count(sym);
  Term prev_base;
  Term state_term;
  if (wire_comb) {
    auto pit = pending_comb_updates_.find(sym);
    if (pit != pending_comb_updates_.end()) {
      prev_base = pit->second;
    } else {
      auto sit = symbol_to_term_.find(sym);
      if (sit != symbol_to_term_.end()) prev_base = sit->second;
    }
  } else if (ctx == StmtContext::NEXT_STATE) {
    auto sit = symbol_to_term_.find(sym);
    if (sit == symbol_to_term_.end()) return;
    state_term = sit->second;
    auto pit = pending_next_updates_.find(state_term);
    prev_base = (pit != pending_next_updates_.end()) ? pit->second : state_term;
  } else {
    // COMBINATIONAL non-wire / INITIAL dynamic-index writes aren't
    // needed by any currently-supported construct; leave unsupported
    // rather than risk an under-constrained partial encoding.
    return;
  }
  if (!prev_base) return;

  Term idx = expr_to_term(sel.selector());
  Term rhs = expr_to_term(rhs_expr);
  rhs = resize_to(rhs, elem_w);
  Term combined = replace_bits_dynamic(prev_base, rhs, idx, elem_w);

  if (wire_comb) {
    if (aliased) pending_comb_aliased_.insert(sym);
    pending_comb_updates_[sym] =
        (condition == solver_->make_term(true))
            ? combined
            : solver_->make_term(Ite, condition, combined, prev_base);
    return;
  }

  // ctx == NEXT_STATE
  pending_next_updates_[state_term] =
      (condition == solver_->make_term(true))
          ? combined
          : solver_->make_term(Ite, condition, combined, prev_base);
}

void SystemVerilogEncoder::refresh_loop_var_term(
    const slang::ast::ValueSymbol & sym)
{
  auto * cur = eval_ctx().findLocal(&sym);
  if (!cur || !cur->isInteger()) {
    throw PonoException("SystemVerilogEncoder: local variable '"
                        + string(sym.name) + "' lost its constant value");
  }
  auto svint = cur->integer();
  uint64_t width = sym.getType().getBitWidth();
  if (width == 0) width = svint.getBitWidth();
  if (width == 0) width = 32;
  Sort sort = solver_->make_sort(BV, width);
  svint.setSigned(false);
  string val_str = svint.toString(slang::LiteralBase::Decimal, false);
  loop_var_terms_[&sym] = solver_->make_term(val_str, sort, 10);
}

void SystemVerilogEncoder::process_statement(const slang::ast::Statement & stmt,
                                             StmtContext ctx,
                                             const Term & condition)
{
  using namespace slang::ast;

  switch (stmt.kind) {
    case StatementKind::ExpressionStatement: {
      auto & es = stmt.as<ExpressionStatement>();
      auto & expr = es.expr;

      // A write to a plain compile-time-unrolled local (a `for`/
      // `while`/`repeat`/`foreach` scratch variable, or any other
      // `VariableDeclaration` local) is neither a wire nor a state
      // variable, so the SMT-term machinery below would silently drop
      // it (falls through to `symbol_to_term_.find() == end()` and
      // returns without doing anything) -- delegate the whole
      // expression to slang's own constant evaluator instead, exactly
      // as `ForLoopStatement`'s own step expressions already do, and
      // refresh the mirrored SMT constant so later condition
      // evaluation (a `while`/`do`-`while` test, a `for`-loop bound)
      // sees the new value.
      {
        const Expression * lval_expr = nullptr;
        if (expr.kind == ExpressionKind::Assignment) {
          lval_expr = &expr.as<AssignmentExpression>().left();
        } else if (expr.kind == ExpressionKind::UnaryOp) {
          auto & unop = expr.as<UnaryExpression>();
          if (unop.op == UnaryOperator::Preincrement
              || unop.op == UnaryOperator::Postincrement
              || unop.op == UnaryOperator::Predecrement
              || unop.op == UnaryOperator::Postdecrement) {
            lval_expr = &unop.operand();
          }
        }
        if (lval_expr) {
          if (auto * base = find_lhs_base(*lval_expr)) {
            auto & vsym = base->as<ValueSymbol>();
            if (eval_ctx().findLocal(&vsym)) {
              if (expr.eval(eval_ctx()).bad()) {
                throw PonoException(
                    "SystemVerilogEncoder: assignment to local variable '"
                    + string(base->name) + "' failed to constant-evaluate");
              }
              refresh_loop_var_term(vsym);
              break;
            }
          }
        }
      }

      // Resolves `lhs_expr` to its base symbol/bit-range/output-alias
      // and the "previous" full-base term a write should compose onto
      // -- shared by plain/compound assignment (whose RHS may read
      // this via an implicit LValueReference) and by `++`/`--` (which
      // reads it directly as "the current value" -- see slice_of()
      // below). Returns nullopt if lhs_expr isn't a shape
      // resolve_lvalue() handles (e.g. a dynamic-index element
      // select), in which case the caller may fall back to
      // process_dynamic_element_assign() for ElementSelect LHSes.
      struct LValueWrite
      {
        const Symbol * sym;
        bool aliased;
        bool has_range;  // true if `sym` itself is an aliased bus slice
        uint64_t lo, hi, slice_w;
        Term prev_base;  // may be null -- see call sites
        bool wire_comb;
        Term state_term;  // only valid when ctx == NEXT_STATE
      };
      auto begin_write =
          [&](const Expression & lhs_expr) -> std::optional<LValueWrite> {
        auto desc = resolve_lvalue(lhs_expr, eval_ctx());
        if (!desc) return std::nullopt;
        const Symbol * sym = desc->base;
        auto resolved = resolve_output_alias(sym);
        sym = resolved.sym;

        uint64_t lo = desc->lo;
        uint64_t hi = desc->hi;
        if (resolved.has_range) {
          // `sym` (pre-alias) is exactly bits [resolved.lo,
          // resolved.hi] of the resolved root -- remap this write
          // into the root's bits.
          lo += resolved.lo;
          hi += resolved.lo;
        }

        LValueWrite w{ sym,
                       resolved.aliased,
                       resolved.has_range,
                       lo,
                       hi,
                       hi - lo + 1,
                       Term(),
                       false,
                       Term() };
        w.wire_comb =
            ctx == StmtContext::COMBINATIONAL && wire_symbols_.count(sym);
        if (w.wire_comb) {
          auto pit = pending_comb_updates_.find(sym);
          if (pit != pending_comb_updates_.end()) w.prev_base = pit->second;
        } else if (ctx == StmtContext::NEXT_STATE) {
          auto sit = symbol_to_term_.find(sym);
          if (sit == symbol_to_term_.end()) return std::nullopt;
          w.state_term = sit->second;
          auto pit = pending_next_updates_.find(w.state_term);
          w.prev_base =
              (pit != pending_next_updates_.end()) ? pit->second : w.state_term;
        } else {
          // COMBINATIONAL non-wire or INITIAL: prev_base is the
          // current (constant) value of the LHS used only for
          // self-reference (compound assignment, ++/--).
          auto sit = symbol_to_term_.find(sym);
          if (sit != symbol_to_term_.end()) w.prev_base = sit->second;
        }
        return w;
      };

      // Slices [lo, hi] out of `base`, or returns a null Term if
      // `base` itself is null (no previous write to reference yet).
      auto slice_of = [&](const Term & base, uint64_t lo, uint64_t hi) -> Term {
        if (!base) return Term();
        uint64_t pw = base->get_sort()->get_width();
        if (lo == 0 && hi == pw - 1) return base;
        return solver_->make_term(Op(Extract, hi, lo), base);
      };

      // Commits `rhs` (already the correct slice_w-wide final value)
      // as the write described by `w` -- shared by plain/compound
      // assignment and by ++/--, which differ only in how `rhs` was
      // computed.
      auto commit_write = [&](const LValueWrite & w, const Term & rhs) {
        if (w.wire_comb) {
          if (w.aliased) pending_comb_aliased_.insert(w.sym);
          // Compose new full-base value from prev_base + slice rhs.
          // On the very first write to this wire within the block
          // there is no prev_base yet; treat it as covering the
          // whole symbol only if it actually does, checked against
          // the symbol's declared width (not the write's own slice
          // width, which is trivially equal to itself) -- otherwise
          // seed a fresh placeholder to splice into, e.g. a
          // `for (i) arr[i] = ...;` pattern writing one element per
          // iteration.
          Term combined;
          if (w.prev_base) {
            combined = replace_bits(w.prev_base, rhs, w.lo, w.hi);
          } else {
            uint64_t sym_w = w.sym->as<ValueSymbol>().getType().getBitWidth();
            bool full_write = !w.has_range && w.lo == 0 && w.hi + 1 == sym_w;
            combined = full_write ? rhs
                                  : replace_bits(
                                        wire_seed_term(w.sym), rhs, w.lo, w.hi);
          }
          if (condition == solver_->make_term(true)) {
            pending_comb_updates_[w.sym] = combined;
          } else {
            Term def = w.prev_base ? w.prev_base : combined;
            pending_comb_updates_[w.sym] =
                solver_->make_term(Ite, condition, combined, def);
          }
          return;
        }

        auto it = symbol_to_term_.find(w.sym);
        if (it == symbol_to_term_.end()) return;
        Term lhs_term = it->second;
        uint64_t base_w = lhs_term->get_sort()->get_width();
        bool full_write = (w.lo == 0 && w.hi == base_w - 1);

        switch (ctx) {
          case StmtContext::NEXT_STATE: {
            Term combined =
                full_write ? rhs : replace_bits(w.prev_base, rhs, w.lo, w.hi);
            Term update;
            if (condition == solver_->make_term(true)) {
              update = combined;
            } else {
              update =
                  solver_->make_term(Ite, condition, combined, w.prev_base);
            }
            pending_next_updates_[w.state_term] = update;
            break;
          }
          case StmtContext::COMBINATIONAL: {
            // Non-wire LHS (e.g. output port reg, or a partially-
            // written base): constrain the appropriate slice via
            // add_constraint, which accepts terms involving input
            // vars (so RHSes that reference input ports work).
            Term lhs_slice =
                full_write
                    ? lhs_term
                    : solver_->make_term(Op(Extract, w.hi, w.lo), lhs_term);
            Term eq = solver_->make_term(Equal, lhs_slice, rhs);
            if (condition != solver_->make_term(true)) {
              eq = solver_->make_term(Implies, condition, eq);
            }
            fts_.add_constraint(eq);
            break;
          }
          case StmtContext::INITIAL: {
            Term lhs_slice =
                full_write
                    ? lhs_term
                    : solver_->make_term(Op(Extract, w.hi, w.lo), lhs_term);
            Term eq = solver_->make_term(Equal, lhs_slice, rhs);
            fts_.constrain_init(eq);
            break;
          }
        }
      };

      if (expr.kind == ExpressionKind::Assignment) {
        auto & assign = expr.as<AssignmentExpression>();
        auto & lhs_expr = assign.left();
        auto & rhs_expr = assign.right();

        auto w = begin_write(lhs_expr);
        if (!w) {
          // resolve_lvalue() only handles constant-index selects; a
          // runtime-variable index (`arr[idx] = rhs`) needs a
          // dynamic-position splice instead of a static bit range.
          if (lhs_expr.kind == ExpressionKind::ElementSelect) {
            process_dynamic_element_assign(
                lhs_expr.as<ElementSelectExpression>(),
                rhs_expr,
                ctx,
                condition);
          }
          break;
        }

        // Stash the slice value for any LValueReference inside rhs.
        Term saved_lvalue = current_lvalue_term_;
        current_lvalue_term_ = slice_of(w->prev_base, w->lo, w->hi);

        Term rhs = expr_to_term(rhs_expr);
        current_lvalue_term_ = saved_lvalue;

        rhs = resize_to(rhs, w->slice_w);
        commit_write(*w, rhs);
      } else if (expr.kind == ExpressionKind::UnaryOp) {
        // `i++`/`--i`/etc. as a standalone statement (distinct from
        // the same operators used as a `for`-loop step expression,
        // which slang's own constant evaluator already handles
        // separately via ForLoopStatement's step evaluation). Per the
        // LRM these are equivalent to `i = i +/- 1`; reuse the exact
        // same lvalue-resolution/commit machinery as plain assignment
        // above, reading the current value directly rather than
        // evaluating an RHS expression.
        auto & unop = expr.as<UnaryExpression>();
        if (unop.op == UnaryOperator::Preincrement
            || unop.op == UnaryOperator::Postincrement
            || unop.op == UnaryOperator::Predecrement
            || unop.op == UnaryOperator::Postdecrement) {
          auto w = begin_write(unop.operand());
          if (!w) {
            logger.log(1,
                       "SystemVerilogEncoder: skipping unsupported ++/-- "
                       "operand shape");
            break;
          }
          Term cur = slice_of(w->prev_base, w->lo, w->hi);
          if (!cur) {
            throw PonoException(
                "SystemVerilogEncoder: '++'/'--' has no previous value to "
                "read for '"
                + std::string(w->sym->name) + "'");
          }
          bool is_inc = unop.op == UnaryOperator::Preincrement
                        || unop.op == UnaryOperator::Postincrement;
          Term one = solver_->make_term(1, cur->get_sort());
          Term new_val = solver_->make_term(is_inc ? BVAdd : BVSub, cur, one);
          commit_write(*w, new_val);
        }
      }
      break;
    }

    case StatementKind::Block: {
      auto & block = stmt.as<BlockStatement>();
      try {
        for_each_stmt_in_block(block, [&](const Statement & s) {
          process_statement(s, ctx, condition);
        });
      }
      catch (const LoopControlSignal & sig) {
        // Absorb a `disable <this block's name>;` reached from inside
        // (stopping the rest of this block); anything else -- a
        // Break/Continue meant for an enclosing ForLoop, or a Disable
        // targeting a different (typically outer) named block --
        // keeps propagating.
        if (sig.kind == LoopControlSignal::Disable && block.blockSymbol
            && sig.disable_target == block.blockSymbol) {
          break;
        }
        throw;
      }
      break;
    }

    case StatementKind::Conditional: {
      auto & cond_stmt = stmt.as<ConditionalStatement>();

      // If the condition is a compile-time constant (e.g. it only
      // references already-unrolled `for`-loop counters), branch on
      // it directly in C++ instead of building a symbolic guard for
      // both arms -- this is what lets break/continue/disable, which
      // can only be modeled as C++-level control flow
      // (LoopControlSignal), propagate correctly out of whichever
      // branch is actually taken. Skipped for the (rare) pattern-match
      // `if` form (multiple conditions); falls back to the general
      // symbolic-guard path below for any condition that isn't
      // const-evaluable (e.g. depends on a runtime signal).
      if (cond_stmt.conditions.size() == 1) {
        auto const_cv = cond_stmt.conditions[0].expr->eval(eval_ctx());
        if (!const_cv.bad()) {
          if (const_cv.isTrue()) {
            process_statement(cond_stmt.ifTrue, ctx, condition);
          } else if (cond_stmt.ifFalse) {
            process_statement(*cond_stmt.ifFalse, ctx, condition);
          }
          break;
        }
      }

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

      // For casex/casez, a constant item pattern's X (casex only) or
      // Z (both; `?` is just an alias for `z` here) bits are
      // wildcards: build a (mask, value) pair with a 0 at each
      // wildcard bit position and 1s everywhere else, so
      // (sel & mask) == value ignores exactly those positions.
      // Returns nullopt for a non-constant pattern (nothing to mask
      // against) or a width too large for this uint64_t-based mask,
      // in which case the caller falls back to plain equality.
      auto casex_mask = [&](const Expression & pat_expr)
          -> std::optional<std::pair<uint64_t, uint64_t>> {
        auto cv = pat_expr.eval(eval_ctx());
        if (!cv.isInteger()) return std::nullopt;
        auto & sv = cv.integer();
        uint64_t w = pat_expr.type->getBitWidth();
        if (w == 0 || w > 64) return std::nullopt;
        uint64_t mask = 0, value = 0;
        for (uint64_t i = 0; i < w; ++i) {
          slang::logic_t bit = sv[static_cast<int32_t>(i)];
          bool wildcard =
              (case_stmt.condition == CaseStatementCondition::WildcardXOrZ)
                  ? bit.isUnknown()
                  : bit.value == slang::logic_t::z.value;
          if (!wildcard) {
            mask |= (uint64_t{ 1 } << i);
            if (bit.value == 1) value |= (uint64_t{ 1 } << i);
          }
        }
        return std::make_pair(mask, value);
      };
      bool is_wildcard_case =
          case_stmt.condition == CaseStatementCondition::WildcardXOrZ
          || case_stmt.condition == CaseStatementCondition::WildcardJustZ;

      Term any_item_matched;
      for (auto & item : case_stmt.items) {
        // Build OR of all patterns matching this item.
        Term item_cond;
        for (auto expr : item.expressions) {
          Term match;
          auto mv = is_wildcard_case ? casex_mask(*expr) : std::nullopt;
          if (mv) {
            uint64_t pat_w = expr->type->getBitWidth();
            Sort pat_sort = solver_->make_sort(BV, pat_w);
            Term mask_term = resize_to(
                solver_->make_term(std::to_string(mv->first), pat_sort, 10),
                sel->get_sort()->get_width());
            Term value_term = resize_to(
                solver_->make_term(std::to_string(mv->second), pat_sort, 10),
                sel->get_sort()->get_width());
            match = solver_->make_term(
                Equal, solver_->make_term(BVAnd, sel, mask_term), value_term);
          } else {
            Term pat = expr_to_term(*expr);
            pat = resize_to(pat, sel->get_sort()->get_width());
            match = solver_->make_term(Equal, sel, pat);
          }
          item_cond =
              item_cond ? solver_->make_term(Or, item_cond, match) : match;
        }
        any_item_matched =
            any_item_matched
                ? solver_->make_term(Or, any_item_matched, item_cond)
                : item_cond;
        Term full_cond = (condition == solver_->make_term(true))
                             ? item_cond
                             : solver_->make_term(And, condition, item_cond);
        process_statement(*item.stmt, ctx, full_cond);
      }
      if (case_stmt.defaultCase) {
        // Default: only when none of the other items matched.
        Term not_matched = any_item_matched
                               ? solver_->make_term(Not, any_item_matched)
                               : solver_->make_term(true);
        Term default_cond =
            (condition == solver_->make_term(true))
                ? not_matched
                : solver_->make_term(And, condition, not_matched);
        process_statement(*case_stmt.defaultCase, ctx, default_cond);
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
      // Handle 'assert', 'assume', 'restrict', and 'cover'.  'assume'/
      // 'restrict' share the exact same property-shape handling as
      // 'assert' below; they differ only in what happens to the
      // resulting boolean once it's built (see the two branches further
      // down). 'cover' is handled via reachability duality (see the
      // design-decision note at the top of this file): `cover
      // property(P)` is checked exactly like `assert property(!P)`, so
      // a "violation" of that surrogate assertion is precisely "P was
      // reached" -- it shares the same safety fast path as assert/
      // assume, just with the boolean negated before the shared
      // disable-window/push logic runs.
      bool is_assumption = ca.assertionKind == AssertionKind::Assume
                           || ca.assertionKind == AssertionKind::Restrict;
      bool is_cover = ca.assertionKind == AssertionKind::CoverProperty;
      if (ca.assertionKind == AssertionKind::Assert || is_assumption
          || is_cover) {
        // Strip the clocking wrapper (the clock event is already
        // baked into our per-cycle abstraction) and any explicit
        // `disable iff` wrapper, recording its condition.  If the
        // statement has no explicit `disable iff`, fall back to the
        // enclosing module's `default disable iff`, if any.
        const AssertionExpr * a = &ca.propertySpec;
        const Expression * disable_expr = nullptr;
        while (true) {
          if (a->kind == AssertionExprKind::Clocking) {
            a = &a->as<ClockingAssertionExpr>().expr;
          } else if (a->kind == AssertionExprKind::DisableIff) {
            auto & di = a->as<DisableIffAssertionExpr>();
            disable_expr = &di.condition;
            a = &di.expr;
          } else {
            break;
          }
        }
        if (!disable_expr && current_scope_) {
          disable_expr = compilation_->getDefaultDisable(*current_scope_);
        }

        Term saved_disable_cond = current_disable_cond_;
        if (disable_expr) {
          Term dc = expr_to_term(*disable_expr);
          current_disable_cond_ = solver_->make_term(
              Distinct, dc, solver_->make_term(0, dc->get_sort()));
        } else {
          current_disable_cond_ = Term();
        }

        // Prefer the pure-safety encoding when the property reduces to
        // a single current-cycle Boolean (plain `assert P`, `always
        // P`, bounded `|->` / `|=>` / `##k` implications).
        // assertion_expr_to_bool returns null as soon as a genuine
        // liveness operator (eventually / unbounded until) appears.
        if (Term prop = assertion_expr_to_bool(*a)) {
          if (is_cover) {
            // Reachability duality: negate before the disable-window
            // exemption below runs, so a `disable iff C` on a cover
            // property correctly means "don't count P as covered while
            // C holds" (the composed surrogate is `!P || C`, whose own
            // violation is exactly "P held and C did not").
            prop = solver_->make_term(Not, prop);
          }
          // Catch-all `disable iff` exemption for property shapes
          // (plain `assert P`, `always P`, ...) that don't already
          // gate themselves more precisely inside assertion_expr_to_bool.
          // For an assumption, this is exactly the right shape too:
          // "assume P disable iff C" means P is only assumed while C is
          // false, i.e. the ever-true constraint is (C || P).
          if (Term dw = disable_window(0)) {
            prop = solver_->make_term(Or, dw, prop);
          }
          if (is_assumption) {
            // Hold at every reachable step (init and, via the transition
            // relation, every subsequent state) -- the same "always true"
            // primitive already used for plain state/input invariants
            // elsewhere in the encoder, just applied to an assumption
            // instead of a proof obligation.
            fts_.add_constraint(prop, /*to_init_and_next=*/true);
            logger.log(1,
                       "SystemVerilogEncoder: extracted assumption "
                       "constraint from {}",
                       make_name(assertion_label(stmt)));
          } else {
            propvec_.push_back(prop);
            logger.log(1,
                       "SystemVerilogEncoder: extracted safety assertion "
                       "property {} (index {})",
                       make_name(assertion_label(stmt)),
                       propvec_.size() - 1);
          }
          current_disable_cond_ = saved_disable_cond;
          break;
        }

        if (is_assumption) {
          // Temporal (non-safety) assume/restrict properties would need
          // their own fairness-constraint machinery (assuming a GF
          // condition rather than proving one), which nothing else in
          // the encoder builds yet -- skip cleanly rather than attempt a
          // partial, likely-unsound translation.
          logger.log(1,
                     "SystemVerilogEncoder: skipping unsupported temporal "
                     "assume/restrict property {}",
                     make_name(assertion_label(stmt)));
          current_disable_cond_ = saved_disable_cond;
          break;
        }

        if (is_cover) {
          // A temporal/sequence-shaped cover goal (e.g. `cover property
          // (a ##1 b)`) would need the reachability duality above
          // extended through the same LTL tableau `ltl_to_sat()` builds
          // for `assert`/`assume` -- negating a liveness obligation
          // doesn't correspond to "was this ever reached" the way it
          // does for a plain current-cycle Boolean, so that extension
          // is out of scope here. Throw rather than silently drop.
          throw PonoException(
              "SystemVerilogEncoder: temporal/sequence-shaped 'cover "
              "property' is not supported");
        }

        // Otherwise build the general LTL tableau for the *negated*
        // property and collect its eventuality-discharge justice
        // conditions.  A fair lasso of the resulting system (every
        // justice condition true infinitely often) on which the
        // negated property holds at cycle 0 is exactly a
        // counterexample to the original assertion.
        TermVec justice;
        Term satpsi = ltl_to_sat(*a, /*neg=*/true, justice);
        if (!satpsi) {
          logger.log(1,
                     "SystemVerilogEncoder: skipping unsupported temporal "
                     "assertion kind {}",
                     static_cast<int>(a->kind));
          current_disable_cond_ = saved_disable_cond;
          break;
        }

        Sort bv1 = solver_->make_sort(BV, 1);
        Term one_bv1 = solver_->make_term(1, bv1);

        // Per-property activation latch: a free 1-bit constant.  The
        // justice set forces it to 1 (so this property's time-0
        // obligation is enabled), while every *other* property's latch
        // may stay 0, leaving their obligations vacuous.  This keeps
        // independent LTL properties from interfering in one system.
        Term act = fts_.make_statevar(
            make_name("__ltl_act_" + std::to_string(latch_counter_++)), bv1);
        fts_.assign_next(act, act);
        Term act_bool = solver_->make_term(Equal, act, one_bv1);

        // Time-0 obligation: when active, the negated property must
        // hold at the first cycle.  Gated by the shared init flag so
        // it constrains only cycle 0, and added to the transition
        // relation (it references the tableau's promise inputs) rather
        // than to the initial-state predicate.
        Term obligation = solver_->make_term(
            Implies,
            solver_->make_term(And, ltl_init_flag(), act_bool),
            satpsi);
        fts_.add_constraint(obligation, /*to_init_and_next=*/false);

        justice.push_back(act_bool);
        ltl_justice_.push_back(justice);
        logger.log(1,
                   "SystemVerilogEncoder: extracted LTL liveness property "
                   "{} (index {}, {} justice condition(s))",
                   make_name(assertion_label(stmt)),
                   ltl_justice_.size() - 1,
                   justice.size());
        current_disable_cond_ = saved_disable_cond;
      }
      break;
    }

    case StatementKind::ImmediateAssertion: {
      // A *procedural* immediate assertion (`assert (expr);`),
      // distinct from the concurrent `assert property (...)` form
      // handled above. Reuses the same Assert-vs-Assume/Restrict
      // split: an assert becomes a safety property, an assume/
      // restrict becomes a standing constraint. Guarded by the
      // accumulated path `condition` (e.g. an enclosing `if`) rather
      // than treated as always-active, since it's only actually
      // reached when program flow gets there -- "if reached, expr
      // must hold" for assert, "if reached, assume expr" for
      // assume/restrict. Pass/fail action blocks (`assert (x) else
      // $error(...);`) are simulation-only display statements with no
      // synthesis meaning and are intentionally not processed, same
      // as $display/$error elsewhere in this encoder. `cover` uses the
      // same reachability-duality contract as the ConcurrentAssertion
      // case above (see the design-decision note at the top of this
      // file): `cover (expr);` is checked exactly like
      // `assert (!expr);`.
      auto & ia = stmt.as<ImmediateAssertionStatement>();
      bool is_assumption = ia.assertionKind == AssertionKind::Assume
                           || ia.assertionKind == AssertionKind::Restrict;
      bool is_cover = ia.assertionKind == AssertionKind::CoverProperty;
      if (ia.assertionKind == AssertionKind::Assert || is_assumption
          || is_cover) {
        Term cond_term = expr_to_term(ia.cond);
        Term bool_cond = solver_->make_term(
            Distinct, cond_term, solver_->make_term(0, cond_term->get_sort()));
        if (is_cover) {
          bool_cond = solver_->make_term(Not, bool_cond);
        }
        Term prop = (condition == solver_->make_term(true))
                        ? bool_cond
                        : solver_->make_term(Implies, condition, bool_cond);
        if (is_assumption) {
          fts_.add_constraint(prop, /*to_init_and_next=*/true);
          logger.log(1,
                     "SystemVerilogEncoder: extracted assumption constraint");
        } else {
          propvec_.push_back(prop);
          logger.log(1,
                     "SystemVerilogEncoder: extracted safety assertion "
                     "property from immediate assertion (index {})",
                     propvec_.size() - 1);
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
        bool broke = false;
        try {
          process_statement(loop.body, ctx, condition);
        }
        catch (const LoopControlSignal & sig) {
          if (sig.kind == LoopControlSignal::Break) {
            broke = true;
          } else if (sig.kind != LoopControlSignal::Continue) {
            // A Disable targeting some other (typically outer) named
            // block keeps propagating past this loop.
            throw;
          }
          // Continue: swallow and fall through to run the step
          // expressions below, matching SV's `continue` (which still
          // runs the step before the next iteration test) -- same as
          // a normally-completed iteration.
        }
        if (broke) break;
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

    case StatementKind::ForeverLoop: {
      // `initial forever @(...) body` (a legacy structural spelling of
      // `always @(...) body`) is recognized and redirected before
      // ever reaching process_statement() at all -- see
      // as_forever_event_body() in process_assignments()/
      // process_instance(). Any `forever` reached *here* is some other
      // shape (no event control, nested inside another statement,
      // etc.): it has no static iteration bound at all, unlike
      // `for`/`while`/`repeat`, which this encoder unrolls up to a
      // compile-time-computable count -- a genuine architectural
      // boundary, not a "not implemented yet" gap. Throw a clear error
      // rather than silently dropping whatever is inside it.
      throw PonoException(
          "SystemVerilogEncoder: 'forever' is only supported as "
          "'initial forever @(...) ...', a structural spelling of "
          "'always @(...) ...'; a bare forever loop has no static "
          "iteration bound and is not supported");
    }

    case StatementKind::WhileLoop: {
      // Compile-time unroll, same contract as `for`: the condition
      // must be constant-evaluable on every iteration (it can
      // reference `for`-loop counters or other locals kept in sync by
      // the ExpressionStatement fast path above) -- a condition that
      // genuinely depends on a runtime (free/registered) signal can't
      // be modeled as C++-level control flow at all, the same
      // architectural boundary as the runtime-dependent break/
      // continue/disable case below.
      auto & loop = stmt.as<WhileLoopStatement>();
      constexpr size_t MAX_ITERS = 65536;
      for (size_t it = 0;; ++it) {
        auto cv = loop.cond.eval(eval_ctx());
        if (cv.bad()) {
          throw PonoException(
              "SystemVerilogEncoder: 'while' condition is not a "
              "compile-time constant (runtime-dependent while loops "
              "are not supported)");
        }
        if (!cv.isTrue()) break;
        if (it >= MAX_ITERS) {
          throw PonoException("SystemVerilogEncoder: 'while' loop exceeded "
                              + std::to_string(MAX_ITERS) + " iterations");
        }
        bool broke = false;
        try {
          process_statement(loop.body, ctx, condition);
        }
        catch (const LoopControlSignal & sig) {
          if (sig.kind == LoopControlSignal::Break) {
            broke = true;
          } else if (sig.kind != LoopControlSignal::Continue) {
            throw;
          }
        }
        if (broke) break;
      }
      break;
    }

    case StatementKind::DoWhileLoop: {
      // Same as WhileLoop, but the condition is tested after the
      // first execution of the body (`do ... while (cond);`).
      auto & loop = stmt.as<DoWhileLoopStatement>();
      constexpr size_t MAX_ITERS = 65536;
      for (size_t it = 0;; ++it) {
        if (it >= MAX_ITERS) {
          throw PonoException("SystemVerilogEncoder: 'do-while' loop exceeded "
                              + std::to_string(MAX_ITERS) + " iterations");
        }
        bool broke = false;
        try {
          process_statement(loop.body, ctx, condition);
        }
        catch (const LoopControlSignal & sig) {
          if (sig.kind == LoopControlSignal::Break) {
            broke = true;
          } else if (sig.kind != LoopControlSignal::Continue) {
            throw;
          }
        }
        if (broke) break;
        auto cv = loop.cond.eval(eval_ctx());
        if (cv.bad()) {
          throw PonoException(
              "SystemVerilogEncoder: 'do-while' condition is not a "
              "compile-time constant (runtime-dependent do-while "
              "loops are not supported)");
        }
        if (!cv.isTrue()) break;
      }
      break;
    }

    case StatementKind::RepeatLoop: {
      // The trip count is evaluated once, up front (per the LRM,
      // `repeat` takes a plain expression, not a re-tested condition
      // like `while`); an unroll-time-unresolvable count (a runtime
      // signal) is out of scope, same contract as `for`/`while`
      // bounds.
      auto & loop = stmt.as<RepeatLoopStatement>();
      auto cv = loop.count.eval(eval_ctx());
      auto n_opt = cv.bad() ? std::nullopt : cv.integer().as<uint64_t>();
      if (!n_opt) {
        throw PonoException(
            "SystemVerilogEncoder: 'repeat' count is not a "
            "compile-time constant (runtime-dependent repeat counts "
            "are not supported)");
      }
      constexpr uint64_t MAX_ITERS = 65536;
      uint64_t n = *n_opt;
      if (n > MAX_ITERS) {
        throw PonoException("SystemVerilogEncoder: 'repeat' loop exceeded "
                            + std::to_string(MAX_ITERS) + " iterations");
      }
      for (uint64_t it = 0; it < n; ++it) {
        bool broke = false;
        try {
          process_statement(loop.body, ctx, condition);
        }
        catch (const LoopControlSignal & sig) {
          if (sig.kind == LoopControlSignal::Break) {
            broke = true;
          } else if (sig.kind != LoopControlSignal::Continue) {
            throw;
          }
        }
        if (broke) break;
      }
      break;
    }

    case StatementKind::ForeachLoop: {
      // Scoped to a single iterated dimension with a concrete loop
      // variable and a statically-known range -- exactly the shape
      // `foreach (arr[i])` produces for a fixed-size packed
      // array/vector, which is all a compile-time-unrolling model can
      // support. Multiple dimensions (`foreach (arr[i][j])`) or a
      // dynamically-sized dimension (a real dynamic array/queue,
      // which this encoder has no sort for anyway) throw a clear
      // error instead of silently iterating the wrong thing or just
      // the first dimension.
      auto & loop = stmt.as<ForeachLoopStatement>();
      if (loop.loopDims.size() != 1 || !loop.loopDims[0].loopVar
          || !loop.loopDims[0].range) {
        throw PonoException(
            "SystemVerilogEncoder: 'foreach' is only supported over a "
            "single statically-sized dimension with a loop variable "
            "(multi-dimensional or dynamically-sized foreach is not "
            "supported)");
      }
      auto & dim = loop.loopDims[0];
      auto & iter_sym = *dim.loopVar;
      int32_t lo = dim.range->lower();
      int32_t hi = dim.range->upper();
      uint64_t width = iter_sym.getType().getBitWidth();
      if (width == 0) width = 32;

      constexpr uint64_t MAX_ITERS = 65536;
      if (static_cast<uint64_t>(hi) - static_cast<uint64_t>(lo) + 1
          > MAX_ITERS) {
        throw PonoException("SystemVerilogEncoder: 'foreach' loop exceeded "
                            + std::to_string(MAX_ITERS) + " iterations");
      }
      for (int32_t idx = lo; idx <= hi; ++idx) {
        slang::SVInt iv(static_cast<slang::bitwidth_t>(width),
                        static_cast<uint64_t>(idx),
                        /*isSigned=*/true);
        eval_ctx().createLocal(&iter_sym, slang::ConstantValue(iv));
        refresh_loop_var_term(iter_sym);
        bool broke = false;
        try {
          process_statement(loop.body, ctx, condition);
        }
        catch (const LoopControlSignal & sig) {
          if (sig.kind != LoopControlSignal::Break
              && sig.kind != LoopControlSignal::Continue) {
            loop_var_terms_.erase(&iter_sym);
            eval_ctx().deleteLocal(&iter_sym);
            throw;
          }
          broke = sig.kind == LoopControlSignal::Break;
        }
        loop_var_terms_.erase(&iter_sym);
        eval_ctx().deleteLocal(&iter_sym);
        if (broke) break;
      }
      break;
    }

    case StatementKind::Break:
    case StatementKind::Continue:
    case StatementKind::Disable: {
      // `condition` is only ever narrowed away from the trivial `true`
      // term by the Conditional/Case cases' *general symbolic* guard
      // building -- the Conditional case's constant-fold fast path
      // (see its comment) always passes `condition` through unchanged,
      // and Block/ForLoop never touch it at all. So reaching this
      // statement with `condition` still exactly `true` means every
      // enclosing `if` along the way was compile-time-resolved, i.e.
      // this control-flow statement is genuinely, unconditionally
      // reached at this point in the unrolling -- interpretable via
      // LoopControlSignal, caught by the nearest enclosing ForLoop
      // (Break/Continue) or matching named Block (Disable). Any other
      // value means we came through at least one runtime-dependent
      // `if`, which (since that path processes both arms
      // unconditionally) can't be correctly modeled as C++-level
      // control flow at all -- a clear error beats silently always- or
      // never-triggering regardless of the real condition.
      if (condition != solver_->make_term(true)) {
        throw PonoException(
            "SystemVerilogEncoder: break/continue/disable is only "
            "supported when its controlling condition is a compile-time "
            "constant (e.g. depends only on already-unrolled for-loop "
            "counters)");
      }
      if (stmt.kind == StatementKind::Break) {
        throw LoopControlSignal{ LoopControlSignal::Break };
      }
      if (stmt.kind == StatementKind::Continue) {
        throw LoopControlSignal{ LoopControlSignal::Continue };
      }
      auto & ds = stmt.as<DisableStatement>();
      const Symbol * target = nullptr;
      if (ds.target.kind == ExpressionKind::ArbitrarySymbol) {
        target = ds.target.as<ArbitrarySymbolExpression>().symbol.get();
      }
      throw LoopControlSignal{ LoopControlSignal::Disable, target };
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

      if (binop.op == BinaryOperator::WildcardEquality
          || binop.op == BinaryOperator::WildcardInequality) {
        // Must special-case *before* the generic eager left/right
        // conversion below: a right operand with unknown bits (e.g.
        // `4'b10??`) would otherwise hit the generic (wildcard-
        // unaware) IntegerLiteral case first and crash trying to hand
        // an X-containing decimal string to the solver, the same bug
        // casex/casez's item patterns had. Per the LRM, only the
        // *right* operand's X/Z bits are wildcards, so left is always
        // safe to convert normally; only convert right if we end up
        // needing it for the plain-equality fallback below. Reuses
        // the same (mask, value) technique as casex/casez in
        // process_statement()'s Case handling: build a mask with a 0
        // at each wildcard bit and compare (left & mask) == value,
        // ignoring exactly those positions. Falls back to plain
        // equality if the right operand isn't a constant (nothing to
        // wildcard against -- this encoder's BV model has no way for
        // a non-literal term to hold an unknown bit at all).
        Term left = expr_to_term(binop.left());
        uint64_t result_width = expr.type->getBitWidth();
        Term eq;
        auto rhs_cv = binop.right().eval(eval_ctx());
        uint64_t rhs_w = binop.right().type->getBitWidth();
        if (!rhs_cv.bad() && rhs_cv.isInteger() && rhs_w > 0 && rhs_w <= 64) {
          auto & sv = rhs_cv.integer();
          uint64_t mask = 0, value = 0;
          for (uint64_t i = 0; i < rhs_w; ++i) {
            slang::logic_t bit = sv[static_cast<int32_t>(i)];
            if (!bit.isUnknown()) {
              mask |= (uint64_t{ 1 } << i);
              if (bit.value == 1) value |= (uint64_t{ 1 } << i);
            }
          }
          Sort rhs_sort = solver_->make_sort(BV, rhs_w);
          Term mask_term =
              resize_to(solver_->make_term(std::to_string(mask), rhs_sort, 10),
                        left->get_sort()->get_width());
          Term value_term =
              resize_to(solver_->make_term(std::to_string(value), rhs_sort, 10),
                        left->get_sort()->get_width());
          eq = solver_->make_term(
              Equal, solver_->make_term(BVAnd, left, mask_term), value_term);
        } else {
          Term right = expr_to_term(binop.right());
          right = resize_to(right, left->get_sort()->get_width());
          eq = solver_->make_term(Equal, left, right);
        }
        Sort bv1 = solver_->make_sort(BV, 1);
        bool want_eq = binop.op == BinaryOperator::WildcardEquality;
        Term result =
            solver_->make_term(Ite,
                               eq,
                               solver_->make_term(want_eq ? 1 : 0, bv1),
                               solver_->make_term(want_eq ? 0 : 1, bv1));
        return resize_to(result, result_width);
      }

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
        case BinaryOperator::Power: {
          // Scoped to a compile-time-constant exponent (the
          // overwhelming majority of real synthesizable use, e.g.
          // `x**2`, `2**WIDTH` in a parameter expression) -- unrolled
          // into repeated multiplication at the operands' common
          // width, matching the truncating-wraparound semantics the
          // rest of this encoder already uses for arithmetic (the
          // uniform resize_to(result, result_width) below then
          // applies exactly as it does for every other operator). A
          // non-constant exponent would need real BV exponentiation,
          // which isn't part of the SMT BV theory and isn't worth a
          // barrel-multiplier-style encoding for how rarely it's used
          // with a runtime exponent.
          auto exp_cv = binop.right().eval(eval_ctx());
          if (exp_cv.bad() || !exp_cv.isInteger()) {
            throw PonoException(
                "SystemVerilogEncoder: '**' is only supported with a "
                "compile-time-constant exponent");
          }
          auto exp_opt = exp_cv.integer().as<uint32_t>();
          if (!exp_opt) {
            throw PonoException(
                "SystemVerilogEncoder: '**' exponent out of range");
          }
          result = solver_->make_term(1, left->get_sort());
          for (uint32_t i = 0; i < *exp_opt; ++i) {
            result = solver_->make_term(BVMul, result, left);
          }
          break;
        }
        case BinaryOperator::CaseEquality: {
          // No X/Z representation in this encoder's pure-BV model, so
          // case equality can never actually differ from logical
          // equality.
          Term eq = solver_->make_term(Equal, left, right);
          Sort bv1 = solver_->make_sort(BV, 1);
          result = solver_->make_term(
              Ite, eq, solver_->make_term(1, bv1), solver_->make_term(0, bv1));
          break;
        }
        case BinaryOperator::CaseInequality: {
          Term eq = solver_->make_term(Equal, left, right);
          Sort bv1 = solver_->make_sort(BV, 1);
          result = solver_->make_term(
              Ite, eq, solver_->make_term(0, bv1), solver_->make_term(1, bv1));
          break;
        }
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
        case UnaryOperator::Plus:
          // Unary `+` is a no-op per the LRM.
          result = operand;
          break;
        case UnaryOperator::BitwiseNand: {
          // Reduction NAND: NOT(AND-reduce) -- same all-ones check as
          // BitwiseAnd above, with the Ite branches swapped.
          Sort bv1 = solver_->make_sort(BV, 1);
          Term all_ones = solver_->make_term(
              Equal,
              operand,
              solver_->make_term(string(operand->get_sort()->get_width(), '1'),
                                 operand->get_sort(),
                                 2));
          result = solver_->make_term(Ite,
                                      all_ones,
                                      solver_->make_term(0, bv1),
                                      solver_->make_term(1, bv1));
          break;
        }
        case UnaryOperator::BitwiseNor: {
          // Reduction NOR: NOT(OR-reduce) -- same any-one check as
          // BitwiseOr above, with the Ite branches swapped.
          Sort bv1 = solver_->make_sort(BV, 1);
          Term any_one = solver_->make_term(
              Distinct, operand, solver_->make_term(0, operand->get_sort()));
          result = solver_->make_term(Ite,
                                      any_one,
                                      solver_->make_term(0, bv1),
                                      solver_->make_term(1, bv1));
          break;
        }
        case UnaryOperator::BitwiseXnor: {
          // Reduction XNOR: NOT(XOR-reduce parity) -- same bit-by-bit
          // XOR fold as BitwiseXor above, negated at the end.
          uint64_t w = operand->get_sort()->get_width();
          Term parity = solver_->make_term(Op(Extract, 0, 0), operand);
          for (uint64_t i = 1; i < w; i++) {
            Term bit = solver_->make_term(Op(Extract, i, i), operand);
            parity = solver_->make_term(BVXor, parity, bit);
          }
          result = solver_->make_term(BVNot, parity);
          break;
        }
        default:
          throw PonoException(
              "SystemVerilogEncoder: unsupported unary operator "
              + to_string(static_cast<int>(unop.op)));
      }
      return resize_to(result, result_width);
    }

    case ExpressionKind::Streaming: {
      // Scoped to a single stream with no `with` sub-range clause --
      // covers real synthesizable usage of streaming concatenation
      // (reversing/regrouping one packed value's bits or byte-lanes,
      // e.g. `{<<{a}}`), not the LRM's full generality (multiple
      // streams, `with` ranges, dynamically-sized queues/arrays).
      auto & sc = expr.as<StreamingConcatenationExpression>();
      if (sc.streams().size() != 1 || sc.streams()[0].withExpr) {
        throw PonoException(
            "SystemVerilogEncoder: unsupported streaming concatenation "
            "shape (only a single stream with no `with` range is "
            "supported)");
      }
      Term base = expr_to_term(*sc.streams()[0].operand);
      uint64_t slice = sc.getSliceSize();
      if (slice == 0) {
        // `{>>{x}}` (left-to-right): a single stream is unchanged.
        return base;
      }
      uint64_t w = base->get_sort()->get_width();
      if (w % slice != 0) {
        throw PonoException(
            "SystemVerilogEncoder: streaming concatenation slice size "
            "does not evenly divide operand width");
      }
      // `{<<{x}}` (right-to-left): reassemble slice-wide chunks in
      // reverse order, taken from x's LSB to MSB.
      Term result;
      for (uint64_t lo = 0; lo < w; lo += slice) {
        Term piece = solver_->make_term(Op(Extract, lo + slice - 1, lo), base);
        result = result ? solver_->make_term(Concat, result, piece) : piece;
      }
      return result;
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

    case ExpressionKind::MemberAccess: {
      // Packed-struct field read (`s.field`): extract the field's bit
      // range, using the FieldSymbol's own bitOffset (LSB-relative,
      // per packed-struct layout) and declared width.
      auto & ma = expr.as<MemberAccessExpression>();
      if (ma.member.kind != SymbolKind::Field) {
        throw PonoException(
            "SystemVerilogEncoder: unsupported member access on "
            + std::string(ma.member.name));
      }
      Term base = expr_to_term(ma.value());
      auto & field = ma.member.as<FieldSymbol>();
      uint64_t w = field.getType().getBitWidth();
      uint64_t lo = field.bitOffset;
      uint64_t hi = lo + w - 1;
      return solver_->make_term(Op(Extract, hi, lo), base);
    }

    case ExpressionKind::StructuredAssignmentPattern: {
      // Packed-struct construction (`'{field: value, ...}`): every
      // field must be given a named setter -- concatenate them MSB
      // first (declaration order), matching packed-struct layout.
      auto & pat = expr.as<StructuredAssignmentPatternExpression>();
      auto & canon = expr.type->getCanonicalType();
      if (canon.kind != SymbolKind::PackedStructType) {
        throw PonoException(
            "SystemVerilogEncoder: unsupported assignment pattern target "
            "type");
      }
      auto & st = canon.as<PackedStructType>();
      std::unordered_map<const Symbol *, const Expression *> setters;
      for (auto & ms : pat.memberSetters) {
        setters[&*ms.member] = &*ms.expr;
      }
      std::vector<Term> parts;
      for (auto & m : st.members()) {
        if (m.kind != SymbolKind::Field) continue;
        auto & field = m.as<FieldSymbol>();
        auto sit = setters.find(&m);
        if (sit == setters.end()) {
          throw PonoException(
              "SystemVerilogEncoder: assignment pattern missing field '"
              + std::string(field.name) + "'");
        }
        Term val = resize_to(expr_to_term(*sit->second),
                             field.getType().getBitWidth());
        parts.push_back(val);
      }
      if (parts.empty()) {
        throw PonoException(
            "SystemVerilogEncoder: assignment pattern for empty struct");
      }
      Term result = parts[0];
      for (size_t i = 1; i < parts.size(); ++i) {
        result = solver_->make_term(Concat, result, parts[i]);
      }
      return result;
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

    case ExpressionKind::Call: {
      // `$past`/`$stable`/`$changed`/`$rose`/`$fell` all expand into
      // (or read from) a chain of 1-cycle latch state vars;
      // `$onehot`/`$onehot0` are a plain bit trick; `$isunknown` is a
      // constant given this encoder's pure 2-valued bitvector model.
      // Other calls (user subroutines, system tasks) are not
      // supported.
      auto & call = expr.as<CallExpression>();
      if (call.isSystemCall() && call.getSubroutineName() == "$past") {
        auto args = call.arguments();
        if (args.empty() || !args[0]) {
          throw PonoException(
              "SystemVerilogEncoder: $past with no value argument");
        }
        Term val = expr_to_term(*args[0]);
        uint32_t n = 1;
        if (args.size() >= 2 && args[1]) {
          auto cv = args[1]->eval(eval_ctx());
          if (cv.isInteger()) {
            auto opt = cv.integer().as<uint32_t>();
            if (opt) n = *opt;
          }
        }
        if (n == 0) return val;
        return make_history_chain(val, n);
      }
      if (call.isSystemCall() && call.getSubroutineName() == "$stable") {
        auto args = call.arguments();
        if (args.empty() || !args[0]) {
          throw PonoException(
              "SystemVerilogEncoder: $stable with no value argument");
        }
        Term val = expr_to_term(*args[0]);
        Term eq = solver_->make_term(Equal, val, make_history_chain(val, 1));
        Sort bv1 = solver_->make_sort(BV, 1);
        return solver_->make_term(
            Ite, eq, solver_->make_term(1, bv1), solver_->make_term(0, bv1));
      }
      if (call.isSystemCall() && call.getSubroutineName() == "$changed") {
        // $stable's negation: the value differs from one cycle ago,
        // over its full width (unlike $rose/$fell, which per the LRM
        // only look at bit 0).
        auto args = call.arguments();
        if (args.empty() || !args[0]) {
          throw PonoException(
              "SystemVerilogEncoder: $changed with no value argument");
        }
        Term val = expr_to_term(*args[0]);
        Term neq =
            solver_->make_term(Distinct, val, make_history_chain(val, 1));
        Sort bv1 = solver_->make_sort(BV, 1);
        return solver_->make_term(
            Ite, neq, solver_->make_term(1, bv1), solver_->make_term(0, bv1));
      }
      if (call.isSystemCall()
          && (call.getSubroutineName() == "$rose"
              || call.getSubroutineName() == "$fell")) {
        // Per the LRM, $rose/$fell only look at bit 0 of a
        // multi-bit argument. Reuses the same 1-cycle latch chain
        // $past/$stable already build for "value one cycle ago".
        bool is_rose = call.getSubroutineName() == "$rose";
        auto args = call.arguments();
        if (args.empty() || !args[0]) {
          throw PonoException("SystemVerilogEncoder: "
                              + std::string(call.getSubroutineName())
                              + " with no value argument");
        }
        Term val = expr_to_term(*args[0]);
        Sort bv1 = solver_->make_sort(BV, 1);
        Term bit0 = solver_->make_term(Op(Extract, 0, 0), val);
        Term prev_bit0 = make_history_chain(bit0, 1);
        Term now_val = solver_->make_term(is_rose ? 1 : 0, bv1);
        Term prev_val = solver_->make_term(is_rose ? 0 : 1, bv1);
        Term edge =
            solver_->make_term(And,
                               solver_->make_term(Equal, bit0, now_val),
                               solver_->make_term(Equal, prev_bit0, prev_val));
        return solver_->make_term(
            Ite, edge, solver_->make_term(1, bv1), solver_->make_term(0, bv1));
      }
      if (call.isSystemCall() && call.getSubroutineName() == "$isunknown") {
        // This encoder's SMT model is pure 2-valued bitvectors -- there
        // is no X/Z representation at all -- so nothing is ever unknown.
        return solver_->make_term(0, solver_->make_sort(BV, 1));
      }
      if (call.isSystemCall()
          && (call.getSubroutineName() == "$onehot"
              || call.getSubroutineName() == "$onehot0")) {
        // Standard power-of-two bit trick: (x & (x-1)) == 0 iff x has
        // at most one bit set; $onehot additionally requires x != 0.
        // No popcount adder needed.
        auto args = call.arguments();
        if (args.empty() || !args[0]) {
          throw PonoException("SystemVerilogEncoder: "
                              + std::string(call.getSubroutineName())
                              + " with no value argument");
        }
        Term val = expr_to_term(*args[0]);
        Term one = solver_->make_term(1, val->get_sort());
        Term minus_one = solver_->make_term(BVSub, val, one);
        Term at_most_one =
            solver_->make_term(Equal,
                               solver_->make_term(BVAnd, val, minus_one),
                               solver_->make_term(0, val->get_sort()));
        Term cond = at_most_one;
        if (call.getSubroutineName() == "$onehot") {
          Term nonzero = solver_->make_term(
              Distinct, val, solver_->make_term(0, val->get_sort()));
          cond = solver_->make_term(And, nonzero, at_most_one);
        }
        Sort bv1 = solver_->make_sort(BV, 1);
        return solver_->make_term(
            Ite, cond, solver_->make_term(1, bv1), solver_->make_term(0, bv1));
      }
      throw PonoException("SystemVerilogEncoder: unsupported call to "
                          + std::string(call.getSubroutineName()));
    }

    default:
      throw PonoException("SystemVerilogEncoder: unsupported expression kind "
                          + to_string(static_cast<int>(expr.kind)));
  }
}

// ============================================================================
// SVA assertion-expression conversion
// ============================================================================

Term SystemVerilogEncoder::make_history_chain(const Term & value, uint32_t n)
{
  Sort sort = value->get_sort();
  Term zero = solver_->make_term(0, sort);
  Term link = value;
  for (uint32_t i = 0; i < n; ++i) {
    Term latch = fts_.make_statevar(
        make_name("__sva_past_" + std::to_string(latch_counter_++)), sort);
    fts_.constrain_init(solver_->make_term(Equal, latch, zero));
    fts_.assign_next(latch, link);
    link = latch;
  }
  return link;
}

namespace {
// Detect a SequenceConcat that we can interpret as a constant
// k-cycle delay applied to a single inner assertion expression
// (`##k Q`).  Returns (k, Q*) on success, std::nullopt otherwise.
std::optional<std::pair<uint32_t, const slang::ast::AssertionExpr *>>
match_const_delay_seq(const slang::ast::AssertionExpr & ae)
{
  using namespace slang::ast;
  if (ae.kind != AssertionExprKind::SequenceConcat) return std::nullopt;
  auto & sc = ae.as<SequenceConcatExpr>();
  if (sc.elements.size() != 1) return std::nullopt;
  auto & e = sc.elements[0];
  if (!e.delay.max || *e.delay.max != e.delay.min) return std::nullopt;
  return std::make_pair(e.delay.min, &*e.sequence);
}

// A named `sequence`/`property` declaration referenced by name (e.g.
// `assert property (p_check);`) binds as a SimpleAssertionExpr wrapping
// an AssertionInstanceExpression -- not a plain boolean Expression --
// so routing it through expr_to_term() throws "unsupported expression
// kind". Slang has already expanded the referenced item's body (with
// any arguments substituted) into `body`; for the no-argument,
// non-recursive case that's exactly the AssertionExpr this encoder
// should recurse into instead. Returns nullptr if `e` isn't such a
// reference (the caller should fall back to its normal expr_to_term()
// path). Argumented, local-variable-bearing, or recursive property
// instantiations are a materially harder problem (they need their own
// binding environment) and are out of scope; throws a clear error
// rather than silently mis-evaluating them.
const slang::ast::AssertionExpr * resolve_named_assertion_ref(
    const slang::ast::Expression & e)
{
  using namespace slang::ast;
  if (e.kind != ExpressionKind::AssertionInstance) return nullptr;
  auto & aie = e.as<AssertionInstanceExpression>();
  if (aie.isRecursiveProperty || !aie.arguments.empty()
      || !aie.localVars.empty()) {
    throw PonoException(
        "SystemVerilogEncoder: named sequence/property references with "
        "arguments, local variables, or recursion are not supported");
  }
  return &aie.body;
}
}  // namespace

// ============================================================================
// LTL tableau
// ============================================================================
//
// Properties that are not pure safety are translated with a standard
// symbolic LTL tableau (temporal testers).  Every temporal operator
// introduces a one-step "promise" tester:
//
//   * a free 1-bit input  n  guessing whether the operator's body
//     holds at the *next* cycle (this is the operator's value now);
//   * a 1-bit state  s  with  s' = n  (it remembers the guess), plus
//     the per-cycle consistency constraint  s == body, which forces
//     the guess made one cycle earlier to have been correct.
//
// Greatest-fixpoint operators (G, R) need no fairness.  Least-fixpoint
// operators (F, strong-U) additionally emit a justice (GF) condition
// that discharges the eventuality, ruling out lassos where a promise
// is made forever but never fulfilled.  ltl_to_sat() pushes negation
// to the leaves on the fly, so the testers it builds always match the
// negation-normal form of the (negated) property.

smt::Term SystemVerilogEncoder::ltl_init_flag()
{
  if (ltl_init_flag_) return ltl_init_flag_;
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  Term flag = fts_.make_statevar(make_name("__ltl_init_flag"), bv1);
  fts_.constrain_init(solver_->make_term(Equal, flag, one_bv1));
  fts_.assign_next(flag, zero_bv1);
  ltl_init_flag_ = solver_->make_term(Equal, flag, one_bv1);
  return ltl_init_flag_;
}

smt::Term SystemVerilogEncoder::ltl_before_cycle(uint32_t k)
{
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  // "is cycle 0" pulse, delayed by 1..k-1 steps gives "is cycle i" for
  // each i in [1, k); their disjunction (with the undelayed pulse) is
  // true exactly during cycles [0, k).
  Term pulse = ltl_init_flag();
  Term result = pulse;
  Term link = solver_->make_term(Ite, pulse, one_bv1, zero_bv1);
  for (uint32_t i = 1; i < k; ++i) {
    Term latch = fts_.make_statevar(
        make_name("__before_cycle_" + std::to_string(latch_counter_++)), bv1);
    fts_.constrain_init(solver_->make_term(Equal, latch, zero_bv1));
    fts_.assign_next(latch, link);
    result = solver_->make_term(
        Or, result, solver_->make_term(Equal, latch, one_bv1));
    link = latch;
  }
  return result;
}

smt::Term SystemVerilogEncoder::disable_window(uint32_t window)
{
  if (!current_disable_cond_) return Term();
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  // OR the disable condition together with its 1..window-cycle-delayed
  // versions, so a `disable iff` that was true anywhere in an
  // antecedent's shift window still exempts the whole property, not
  // just the single cycle the check ends up anchored at.
  Term result = current_disable_cond_;
  Term link = solver_->make_term(Ite, current_disable_cond_, one_bv1, zero_bv1);
  for (uint32_t i = 0; i < window; ++i) {
    Term latch = fts_.make_statevar(
        make_name("__disable_hist_" + std::to_string(latch_counter_++)), bv1);
    fts_.constrain_init(solver_->make_term(Equal, latch, zero_bv1));
    fts_.assign_next(latch, link);
    result = solver_->make_term(
        Or, result, solver_->make_term(Equal, latch, one_bv1));
    link = latch;
  }
  return result;
}

smt::Term SystemVerilogEncoder::ltl_make_X(const smt::Term & phi)
{
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  uint32_t id = latch_counter_++;
  Term n = fts_.make_inputvar(make_name("__ltl_n_" + std::to_string(id)), bv1);
  Term s = fts_.make_statevar(make_name("__ltl_s_" + std::to_string(id)), bv1);
  fts_.assign_next(s, n);  // s@(t+1) = n@t  (remember the guess)
  // The guess made at t-1 about phi@t must have been correct.  The
  // cycle-0 instance only pins the otherwise-unused initial s.
  Term phi_bv = solver_->make_term(Ite, phi, one_bv1, zero_bv1);
  fts_.add_constraint(solver_->make_term(Equal, s, phi_bv),
                      /*to_init_and_next=*/false);
  return solver_->make_term(Equal, n, one_bv1);  // sat(X phi) = guess
}

smt::Term SystemVerilogEncoder::ltl_make_G(const smt::Term & phi)
{
  // G phi == phi && X(G phi)
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  uint32_t id = latch_counter_++;
  Term n = fts_.make_inputvar(make_name("__ltl_n_" + std::to_string(id)), bv1);
  Term s = fts_.make_statevar(make_name("__ltl_s_" + std::to_string(id)), bv1);
  fts_.assign_next(s, n);
  Term n_bool = solver_->make_term(Equal, n, one_bv1);
  Term body = solver_->make_term(And, phi, n_bool);  // sat(G phi)
  Term body_bv = solver_->make_term(Ite, body, one_bv1, zero_bv1);
  fts_.add_constraint(solver_->make_term(Equal, s, body_bv),
                      /*to_init_and_next=*/false);
  return body;
}

smt::Term SystemVerilogEncoder::ltl_make_F(const smt::Term & phi,
                                           smt::TermVec & justice)
{
  // F phi == phi || X(F phi)
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  uint32_t id = latch_counter_++;
  Term n = fts_.make_inputvar(make_name("__ltl_n_" + std::to_string(id)), bv1);
  Term s = fts_.make_statevar(make_name("__ltl_s_" + std::to_string(id)), bv1);
  fts_.assign_next(s, n);
  Term n_bool = solver_->make_term(Equal, n, one_bv1);
  Term body = solver_->make_term(Or, phi, n_bool);  // sat(F phi)
  Term body_bv = solver_->make_term(Ite, body, one_bv1, zero_bv1);
  fts_.add_constraint(solver_->make_term(Equal, s, body_bv),
                      /*to_init_and_next=*/false);
  // Discharge: infinitely often phi holds or we stop promising F.
  justice.push_back(
      solver_->make_term(Or, phi, solver_->make_term(Not, n_bool)));
  return body;
}

smt::Term SystemVerilogEncoder::ltl_make_R(const smt::Term & a,
                                           const smt::Term & b)
{
  // a R b == b && (a || X(a R b))   (greatest fixpoint, no fairness)
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  uint32_t id = latch_counter_++;
  Term n = fts_.make_inputvar(make_name("__ltl_n_" + std::to_string(id)), bv1);
  Term s = fts_.make_statevar(make_name("__ltl_s_" + std::to_string(id)), bv1);
  fts_.assign_next(s, n);
  Term n_bool = solver_->make_term(Equal, n, one_bv1);
  Term body = solver_->make_term(
      And, b, solver_->make_term(Or, a, n_bool));  // sat(a R b)
  Term body_bv = solver_->make_term(Ite, body, one_bv1, zero_bv1);
  fts_.add_constraint(solver_->make_term(Equal, s, body_bv),
                      /*to_init_and_next=*/false);
  return body;
}

smt::Term SystemVerilogEncoder::ltl_make_U(const smt::Term & a,
                                           const smt::Term & b,
                                           smt::TermVec & justice)
{
  // a U b (strong) == b || (a && X(a U b))
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  uint32_t id = latch_counter_++;
  Term n = fts_.make_inputvar(make_name("__ltl_n_" + std::to_string(id)), bv1);
  Term s = fts_.make_statevar(make_name("__ltl_s_" + std::to_string(id)), bv1);
  fts_.assign_next(s, n);
  Term n_bool = solver_->make_term(Equal, n, one_bv1);
  Term body = solver_->make_term(
      Or, b, solver_->make_term(And, a, n_bool));  // sat(a U b)
  Term body_bv = solver_->make_term(Ite, body, one_bv1, zero_bv1);
  fts_.add_constraint(solver_->make_term(Equal, s, body_bv),
                      /*to_init_and_next=*/false);
  // Discharge: infinitely often b holds or a U b is false.
  justice.push_back(solver_->make_term(Or, solver_->make_term(Not, body), b));
  return body;
}

Term SystemVerilogEncoder::delay_bool(const Term & cond, uint32_t n)
{
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  Term cond_bv = solver_->make_term(Ite, cond, one_bv1, zero_bv1);
  Term delayed_bv = make_history_chain(cond_bv, n);
  return solver_->make_term(Equal, delayed_bv, one_bv1);
}

namespace {
// The largest total span (in cycles) offsets_ending_now() will build
// before giving up -- a defensive cap against a pathological/absurd
// bounded sequence, mirroring the MAX_ITERS-style caps used elsewhere
// in this encoder for other compile-time-unrolled constructs.
constexpr uint32_t MAX_SEQ_WINDOW = 256;
}  // namespace

smt::TermVec SystemVerilogEncoder::offsets_ending_now(
    const slang::ast::AssertionExpr & seq)
{
  using namespace slang::ast;

  // A single Boolean expression, optionally with a consecutive
  // repetition (`expr[*n:m]`, `expr[+]`, `expr[*]`). Shared by both
  // SimpleAssertionExpr and SequenceWithMatchExpr, which each carry
  // their own std::optional<SequenceRepetition>.
  auto boolean_with_repetition =
      [&](const slang::ast::Expression & expr,
          const std::optional<SequenceRepetition> & repetition) -> TermVec {
    Term b = expr_to_term(expr);
    b = solver_->make_term(Distinct, b, solver_->make_term(0, b->get_sort()));
    if (!repetition) return { b };
    if (repetition->kind != SequenceRepetition::Consecutive) {
      throw PonoException(
          "SystemVerilogEncoder: nonconsecutive/goto sequence repetition "
          "([->]/[=]) is not supported");
    }
    if (!repetition->range.max) {
      throw PonoException(
          "SystemVerilogEncoder: unbounded consecutive sequence "
          "repetition ([*]/[+]/[*n:$]) is not supported");
    }
    uint32_t lo = repetition->range.min;
    uint32_t hi = *repetition->range.max;
    if (hi >= MAX_SEQ_WINDOW) {
      throw PonoException("SystemVerilogEncoder: sequence repetition exceeds "
                          + std::to_string(MAX_SEQ_WINDOW) + " cycles");
    }
    // expr[*lo:hi] matches (ending now, started L cycles ago) for
    // L in [lo-1, hi-1], requiring expr to hold at every one of the
    // L+1 cycles from L-cycles-ago through now.
    TermVec out(hi, Term());
    Term running = b;
    for (uint32_t count = 1; count <= hi; ++count) {
      if (count > 1) {
        running = solver_->make_term(And, running, delay_bool(b, count - 1));
      }
      if (count >= lo) out[count - 1] = running;
    }
    return out;
  };

  // ORs together every non-null entry of a TermVec (the "does it
  // complete here at all" merge, without recomputing offsets_ending_now
  // -- used by Within/Intersect below).
  auto or_vec = [&](const TermVec & v) -> Term {
    Term result;
    for (auto & t : v) {
      if (!t) continue;
      result = result ? solver_->make_term(Or, result, t) : t;
    }
    return result;
  };

  // OR of `base` and its delayed copies over the last `k` cycles
  // (`base` itself, plus 1..k cycles ago) -- "did `base` hold at *some*
  // point in the last k+1 cycles". Used by Within.
  auto window_or = [&](const Term & base, uint32_t k) -> Term {
    Term result = base;
    for (uint32_t j = 1; j <= k; ++j) {
      result = solver_->make_term(Or, result, delay_bool(base, j));
    }
    return result;
  };

  // AND of `base` and its delayed copies over the last `k` cycles --
  // "did `base` hold at *every* point in the last k+1 cycles". Used by
  // Throughout.
  auto window_and = [&](const Term & base, uint32_t k) -> Term {
    Term result = base;
    for (uint32_t j = 1; j <= k; ++j) {
      result = solver_->make_term(And, result, delay_bool(base, j));
    }
    return result;
  };

  switch (seq.kind) {
    case AssertionExprKind::Simple: {
      auto & simple = seq.as<SimpleAssertionExpr>();
      return boolean_with_repetition(simple.expr, simple.repetition);
    }

    case AssertionExprKind::SequenceWithMatch: {
      auto & swm = seq.as<SequenceWithMatchExpr>();
      // A parenthesized sequence with its own repetition
      // (`(seq)[*n:m]`) would need to convolve `seq`'s own offset
      // vector with itself count times -- out of scope; only a plain
      // Boolean operand with repetition is handled today.
      if (swm.repetition) {
        if (swm.expr.kind != AssertionExprKind::Simple
            || swm.expr.as<SimpleAssertionExpr>().repetition) {
          throw PonoException(
              "SystemVerilogEncoder: repetition of a non-Boolean sequence "
              "is not supported");
        }
        return boolean_with_repetition(swm.expr.as<SimpleAssertionExpr>().expr,
                                       swm.repetition);
      }
      return offsets_ending_now(swm.expr);
    }

    case AssertionExprKind::FirstMatch:
      // first_match(seq) only restricts *which* match is reported when
      // a sequence can match in more than one way -- it never changes
      // whether a match exists at all, which is all this encoder's
      // callers (an implication antecedent, an intersect/within/
      // throughout operand) ever ask offsets_ending_now() for.
      return offsets_ending_now(seq.as<FirstMatchAssertionExpr>().seq);

    case AssertionExprKind::Clocking:
      // Per this file's multiclock design decision: every named clock
      // is treated as the same global pono-cycle, so a nested clocking
      // change inside a sequence element is simply unwrapped.
      return offsets_ending_now(seq.as<ClockingAssertionExpr>().expr);

    case AssertionExprKind::SequenceConcat: {
      auto & sc = seq.as<SequenceConcatExpr>();
      TermVec acc;
      for (size_t i = 0; i < sc.elements.size(); ++i) {
        auto & elem = sc.elements[i];
        if (!elem.delay.max) {
          throw PonoException(
              "SystemVerilogEncoder: unbounded sequence delay (##[m:$]) is "
              "not supported");
        }
        uint32_t dmin = elem.delay.min;
        uint32_t dmax = *elem.delay.max;
        TermVec elem_offsets = offsets_ending_now(*elem.sequence);
        if (elem_offsets.empty()) return {};

        if (i == 0) {
          // The delay before the very first element just relabels how
          // far back "the sequence's start" is, with no extra
          // condition to AND in.
          size_t new_size = elem_offsets.size() - 1 + dmax + 1;
          if (new_size > MAX_SEQ_WINDOW) {
            throw PonoException("SystemVerilogEncoder: sequence window exceeds "
                                + std::to_string(MAX_SEQ_WINDOW) + " cycles");
          }
          acc.assign(new_size, Term());
          for (size_t l = 0; l < elem_offsets.size(); ++l) {
            if (!elem_offsets[l]) continue;
            for (uint32_t d = dmin; d <= dmax; ++d) {
              size_t idx = l + d;
              acc[idx] = acc[idx]
                             ? solver_->make_term(Or, acc[idx], elem_offsets[l])
                             : elem_offsets[l];
            }
          }
          continue;
        }

        size_t new_size = acc.size() - 1 + dmax + (elem_offsets.size() - 1) + 1;
        if (new_size > MAX_SEQ_WINDOW) {
          throw PonoException("SystemVerilogEncoder: sequence window exceeds "
                              + std::to_string(MAX_SEQ_WINDOW) + " cycles");
        }
        TermVec new_acc(new_size, Term());
        for (size_t lp = 0; lp < acc.size(); ++lp) {
          if (!acc[lp]) continue;
          for (uint32_t d = dmin; d <= dmax; ++d) {
            for (size_t le = 0; le < elem_offsets.size(); ++le) {
              if (!elem_offsets[le]) continue;
              size_t idx = lp + d + le;
              // Bring the (already-anchored-at-"now") prefix condition
              // back by (d + le) cycles so it aligns with this
              // element's own occurrence, then AND with this
              // element's own (unshifted) completion condition.
              Term shifted_prefix =
                  (d + le == 0) ? acc[lp] : delay_bool(acc[lp], d + le);
              Term combined =
                  solver_->make_term(And, shifted_prefix, elem_offsets[le]);
              new_acc[idx] =
                  new_acc[idx] ? solver_->make_term(Or, new_acc[idx], combined)
                               : combined;
            }
          }
        }
        acc = std::move(new_acc);
      }
      return acc;
    }

    case AssertionExprKind::Binary: {
      auto & b = seq.as<BinaryAssertionExpr>();
      switch (b.op) {
        case BinaryAssertionOperator::Intersect: {
          // s1 intersect s2 matches iff both match with the *same*
          // span -- AND the two offset vectors entry-by-entry.
          TermVec v1 = offsets_ending_now(b.left);
          TermVec v2 = offsets_ending_now(b.right);
          if (v1.empty() || v2.empty()) return {};
          size_t n = std::min(v1.size(), v2.size());
          TermVec out(n, Term());
          for (size_t k = 0; k < n; ++k) {
            if (v1[k] && v2[k]) out[k] = solver_->make_term(And, v1[k], v2[k]);
          }
          return out;
        }

        case BinaryAssertionOperator::Within: {
          // s1 within s2 matches over the same span as a match of s2,
          // provided s1 matches ending somewhere inside that span: for
          // each of s2's own completion offsets k, "s1 matched
          // somewhere in the last k+1 cycles" is exactly window_or()
          // over s1's merged "matches here" term.
          TermVec v1 = offsets_ending_now(b.left);
          TermVec v2 = offsets_ending_now(b.right);
          if (v1.empty() || v2.empty()) return {};
          Term s1_matches = or_vec(v1);
          if (!s1_matches) return {};
          TermVec out(v2.size(), Term());
          for (size_t k = 0; k < v2.size(); ++k) {
            if (!v2[k]) continue;
            out[k] = solver_->make_term(
                And, v2[k], window_or(s1_matches, (uint32_t)k));
          }
          return out;
        }

        case BinaryAssertionOperator::Throughout: {
          // expr throughout seq: the plain boolean expr must hold at
          // every cycle spanned by seq's match -- for each of seq's
          // own completion offsets k, that span is exactly the last
          // k+1 cycles, checked via window_and().
          Term expr_bool = assertion_expr_to_bool(b.left);
          if (!expr_bool) return {};
          TermVec v2 = offsets_ending_now(b.right);
          if (v2.empty()) return {};
          TermVec out(v2.size(), Term());
          for (size_t k = 0; k < v2.size(); ++k) {
            if (!v2[k]) continue;
            out[k] = solver_->make_term(
                And, v2[k], window_and(expr_bool, (uint32_t)k));
          }
          return out;
        }

        default:
          // And/Or/Iff/Implies/Until*/FollowedBy/etc. as a *sequence*
          // operand aren't sequence-composition operators in the SVA
          // sense (they combine boolean/property values, not match
          // spans) -- not something offsets_ending_now() is ever
          // asked for by this encoder's callers today.
          return {};
      }
    }

    default:
      // Unsupported sequence shape (a nested StrongWeak/etc. operand
      // this primitive doesn't model yet) -- the caller falls back to
      // its existing unsupported-construct handling.
      return {};
  }
}

smt::Term SystemVerilogEncoder::match_exists(
    const slang::ast::AssertionExpr & seq)
{
  TermVec offsets = offsets_ending_now(seq);
  Term result;
  for (auto & t : offsets) {
    if (!t) continue;
    result = result ? solver_->make_term(Or, result, t) : t;
  }
  return result;
}

smt::Term SystemVerilogEncoder::leading_condition(
    const slang::ast::AssertionExpr & seq)
{
  using namespace slang::ast;
  switch (seq.kind) {
    case AssertionExprKind::Simple: {
      auto & simple = seq.as<SimpleAssertionExpr>();
      if (simple.repetition) {
        throw PonoException(
            "SystemVerilogEncoder: weak()/strong() of a sequence with its "
            "own leading repetition is not supported");
      }
      Term t = expr_to_term(simple.expr);
      return solver_->make_term(
          Distinct, t, solver_->make_term(0, t->get_sort()));
    }
    case AssertionExprKind::FirstMatch:
      return leading_condition(seq.as<FirstMatchAssertionExpr>().seq);
    case AssertionExprKind::Clocking:
      return leading_condition(seq.as<ClockingAssertionExpr>().expr);
    case AssertionExprKind::SequenceConcat:
      return leading_condition(
          *seq.as<SequenceConcatExpr>().elements[0].sequence);
    default:
      throw PonoException(
          "SystemVerilogEncoder: weak()/strong() of this sequence shape is "
          "not supported");
  }
}

smt::Term SystemVerilogEncoder::weak_seq_bool(
    const slang::ast::AssertionExpr & seq)
{
  TermVec offsets = offsets_ending_now(seq);
  if (offsets.empty()) return Term();
  Term me;
  for (auto & t : offsets) {
    if (!t) continue;
    me = me ? solver_->make_term(Or, me, t) : t;
  }
  if (!me) return Term();

  // S = the sequence's own maximum span: the last possible cycle an
  // attempt that started here could still complete by.
  uint32_t s = static_cast<uint32_t>(offsets.size()) - 1;
  Term started_s_ago = delay_bool(leading_condition(seq), s);
  Term completed_in_window = me;
  for (uint32_t j = 1; j <= s; ++j) {
    completed_in_window =
        solver_->make_term(Or, completed_in_window, delay_bool(me, j));
  }
  // Violated iff an attempt began exactly S cycles ago and no
  // completion happened anywhere from then through now; weak(seq) is
  // the negation -- no obligation to ever attempt, but an attempt that
  // did begin must not be a definite, provable failure.
  Term violated = solver_->make_term(
      And, started_s_ago, solver_->make_term(Not, completed_in_window));
  return solver_->make_term(Not, violated);
}

smt::Term SystemVerilogEncoder::ltl_to_sat(const slang::ast::AssertionExpr & ae,
                                           bool neg,
                                           smt::TermVec & justice)
{
  using namespace slang::ast;

  switch (ae.kind) {
    case AssertionExprKind::Clocking:
      return ltl_to_sat(ae.as<ClockingAssertionExpr>().expr, neg, justice);

    case AssertionExprKind::StrongWeak: {
      auto & sw = ae.as<StrongWeakAssertionExpr>();
      // The strong/weak qualifier only meaningfully differs for a
      // genuine bounded sequence -- must it eventually complete
      // (strong) or not (weak)? Any other shape (already a plain
      // Boolean/temporal expression) is unaffected by the qualifier
      // under this encoder's infinite-lasso semantics; just unwrap.
      if (sw.strength == StrongWeakAssertionExpr::Strong) {
        Term me = match_exists(sw.expr);
        if (me) {
          // strong(seq): a genuine liveness obligation -- the sequence
          // must eventually complete a match.
          return neg ? ltl_make_G(solver_->make_term(Not, me))
                     : ltl_make_F(me, justice);
        }
      }
      return ltl_to_sat(sw.expr, neg, justice);
    }

    case AssertionExprKind::Simple: {
      auto & simple = ae.as<SimpleAssertionExpr>();
      if (auto * named = resolve_named_assertion_ref(simple.expr)) {
        return ltl_to_sat(*named, neg, justice);
      }
      if (simple.repetition) {
        // See the matching check in assertion_expr_to_bool(): route
        // through the general bounded sequence matcher (which throws
        // for an unbounded repeat count) instead of silently ignoring
        // the repetition.
        Term me = match_exists(ae);
        if (!me) return Term();
        return neg ? solver_->make_term(Not, me) : me;
      }
      Term t = expr_to_term(simple.expr);
      if (!t) return Term();
      Term zero = solver_->make_term(0, t->get_sort());
      Term b = solver_->make_term(Distinct, t, zero);
      return neg ? solver_->make_term(Not, b) : b;
    }

    case AssertionExprKind::SequenceConcat: {
      // A bare `##k Q` property has the same truth value as Q under
      // our infinite-time semantics (modulo a front shift).  Unwrap.
      if (auto m = match_const_delay_seq(ae)) {
        return ltl_to_sat(*m->second, neg, justice);
      }
      return Term();
    }

    case AssertionExprKind::Unary: {
      auto & u = ae.as<UnaryAssertionExpr>();
      switch (u.op) {
        case UnaryAssertionOperator::Not:
          return ltl_to_sat(u.expr, !neg, justice);

        case UnaryAssertionOperator::Always:
        case UnaryAssertionOperator::SAlways: {
          // G phi  (positive)  /  !G phi == F !phi  (negated)
          Term phi = ltl_to_sat(u.expr, neg, justice);
          if (!phi) return Term();
          return neg ? ltl_make_F(phi, justice) : ltl_make_G(phi);
        }

        case UnaryAssertionOperator::Eventually:
        case UnaryAssertionOperator::SEventually: {
          // F phi  (positive)  /  !F phi == G !phi  (negated)
          Term phi = ltl_to_sat(u.expr, neg, justice);
          if (!phi) return Term();
          return neg ? ltl_make_G(phi) : ltl_make_F(phi, justice);
        }

        case UnaryAssertionOperator::NextTime:
        case UnaryAssertionOperator::SNextTime: {
          // !X phi == X !phi, so the negation rides along inside phi.
          Term phi = ltl_to_sat(u.expr, neg, justice);
          if (!phi) return Term();
          return ltl_make_X(phi);
        }

        default: return Term();
      }
    }

    case AssertionExprKind::Binary: {
      auto & b = ae.as<BinaryAssertionExpr>();
      switch (b.op) {
        case BinaryAssertionOperator::And:
        case BinaryAssertionOperator::Or: {
          Term l = ltl_to_sat(b.left, neg, justice);
          Term r = ltl_to_sat(b.right, neg, justice);
          if (!l || !r) return Term();
          bool is_and = (b.op == BinaryAssertionOperator::And);
          if (neg) is_and = !is_and;  // De Morgan
          return solver_->make_term(is_and ? And : Or, l, r);
        }

        case BinaryAssertionOperator::Iff: {
          // Children appear in both polarities, so build them
          // positively and apply the outer negation here.  (Temporal
          // operands under `iff` are not negation-normalized; such
          // properties are exotic and out of scope for fairness.)
          Term l = ltl_to_sat(b.left, false, justice);
          Term r = ltl_to_sat(b.right, false, justice);
          if (!l || !r) return Term();
          return neg ? solver_->make_term(Distinct, l, r)
                     : solver_->make_term(Equal, l, r);
        }

        case BinaryAssertionOperator::Implies: {
          // a implies b == !a || b ;  !(a implies b) == a && !b
          Term l = ltl_to_sat(b.left, !neg, justice);
          Term r = ltl_to_sat(b.right, neg, justice);
          if (!l || !r) return Term();
          return solver_->make_term(neg ? And : Or, l, r);
        }

        case BinaryAssertionOperator::OverlappedImplication:
        case BinaryAssertionOperator::NonOverlappedImplication: {
          // seq |-> prop / seq |=> prop with a Boolean antecedent:
          //   !a || X^delay b   (delay 0 overlapped, 1 non-overlapped,
          //                       plus any ##k on the consequent)
          // and its negation a && X^delay !b.
          uint32_t delay =
              (b.op == BinaryAssertionOperator::NonOverlappedImplication) ? 1
                                                                          : 0;
          const AssertionExpr * rhs = &b.right;
          if (auto m = match_const_delay_seq(b.right)) {
            delay += m->first;
            rhs = m->second;
          }
          Term l = ltl_to_sat(b.left, !neg, justice);
          Term r = ltl_to_sat(*rhs, neg, justice);
          if (!l || !r) return Term();
          for (uint32_t i = 0; i < delay; ++i) r = ltl_make_X(r);
          return solver_->make_term(neg ? And : Or, l, r);
        }

        case BinaryAssertionOperator::Until:
        case BinaryAssertionOperator::SUntil:
        case BinaryAssertionOperator::UntilWith:
        case BinaryAssertionOperator::SUntilWith: {
          bool strong = (b.op == BinaryAssertionOperator::SUntil
                         || b.op == BinaryAssertionOperator::SUntilWith);
          bool with = (b.op == BinaryAssertionOperator::UntilWith
                       || b.op == BinaryAssertionOperator::SUntilWith);
          if (!neg) {
            Term l = ltl_to_sat(b.left, false, justice);
            Term r = ltl_to_sat(b.right, false, justice);
            if (!l || !r) return Term();
            // until_with: the terminating cycle must also satisfy a.
            Term term = with ? solver_->make_term(And, l, r) : r;
            if (strong) return ltl_make_U(l, term, justice);
            // weak until: a W term == term R (a || term).
            return ltl_make_R(term, solver_->make_term(Or, l, term));
          }
          // Negated, with operands already in negation-normal form
          // (nl = sat(!left), nr = sat(!right)):
          //   !(a U_strong b)      = !a R !b
          //   !(a W b)             = !b U (!a && !b)
          //   !(a U_strong (a&&b)) = !a R (!a || !b)
          //   !(a W (a&&b))        = (!a || !b) U !a
          Term nl = ltl_to_sat(b.left, true, justice);
          Term nr = ltl_to_sat(b.right, true, justice);
          if (!nl || !nr) return Term();
          if (!with) {
            if (strong) return ltl_make_R(nl, nr);
            return ltl_make_U(nr, solver_->make_term(And, nl, nr), justice);
          }
          Term nterm = solver_->make_term(Or, nl, nr);  // !(a && b)
          if (strong) return ltl_make_R(nl, nterm);
          return ltl_make_U(nterm, nl, justice);
        }

        default:
          // Intersect / Throughout / Within / FollowedBy: multi-cycle
          // sequence operators the tableau does not model.
          return Term();
      }
    }

    default: return Term();
  }
}

Term SystemVerilogEncoder::assertion_expr_to_bool(
    const slang::ast::AssertionExpr & ae)
{
  using namespace slang::ast;

  switch (ae.kind) {
    case AssertionExprKind::Clocking: {
      // The clocking event has already been baked into our cycle
      // abstraction; just recurse into the underlying expression.
      return assertion_expr_to_bool(ae.as<ClockingAssertionExpr>().expr);
    }

    case AssertionExprKind::Simple: {
      auto & simple = ae.as<SimpleAssertionExpr>();
      if (auto * named = resolve_named_assertion_ref(simple.expr)) {
        return assertion_expr_to_bool(*named);
      }
      if (simple.repetition) {
        // `expr[*n:m]`/`expr[+]`/`expr[*]`: route through the general
        // bounded sequence matcher (which throws for an unbounded
        // repeat count) instead of silently ignoring the repetition
        // and returning a bare `bool(expr)`.
        return match_exists(ae);
      }
      Term t = expr_to_term(simple.expr);
      // Normalize to Bool: t != 0.
      Sort sort = t->get_sort();
      Term zero = solver_->make_term(0, sort);
      return solver_->make_term(Distinct, t, zero);
    }

    case AssertionExprKind::SequenceConcat: {
      // A standalone `##k Q` as a property: in our infinite-time
      // safety encoding the constant front-shift doesn't change the
      // truth value (it just postpones when the first violation can
      // be reported), so unwrap the inner sequence.
      if (auto matched = match_const_delay_seq(ae)) {
        return assertion_expr_to_bool(*matched->second);
      }
      return Term();
    }

    case AssertionExprKind::StrongWeak: {
      auto & sw = ae.as<StrongWeakAssertionExpr>();
      if (sw.strength == StrongWeakAssertionExpr::Strong) {
        // strong(seq) is a genuine liveness obligation ("must
        // eventually complete"), not reducible to a current-cycle
        // Boolean -- handled by ltl_to_sat()'s StrongWeak case
        // instead, which builds the eventuality tableau. Returning
        // null here forces the ConcurrentAssertion handler to fall
        // through to that path rather than (incorrectly) treating a
        // strong sequence as if it were always-true right now.
        return Term();
      }
      // weak(seq): no obligation to ever match, but an attempt that
      // did begin must not be a definite, provable failure -- see
      // weak_seq_bool(). Any other shape (already a plain Boolean/
      // temporal expression) is unaffected by the qualifier; just
      // unwrap.
      if (Term w = weak_seq_bool(sw.expr)) return w;
      return assertion_expr_to_bool(sw.expr);
    }

    case AssertionExprKind::Unary: {
      auto & u = ae.as<UnaryAssertionExpr>();
      switch (u.op) {
        case UnaryAssertionOperator::Not: {
          Term inner = assertion_expr_to_bool(u.expr);
          if (!inner) return Term();
          return solver_->make_term(Not, inner);
        }
        case UnaryAssertionOperator::Always:
        case UnaryAssertionOperator::SAlways: {
          // Pure safety: `always P` is true at the current cycle
          // exactly when P is true at the current cycle (the
          // "always at every cycle" closure is implicit in the
          // per-cycle property check).
          return assertion_expr_to_bool(u.expr);
        }
        default:
          // Eventually / SEventually / NextTime / SNextTime can't
          // be folded into a current-cycle Boolean: liveness must
          // come through the top-level dispatch and NextTime needs
          // a forward-shift the encoder doesn't model yet.
          return Term();
      }
    }

    case AssertionExprKind::Binary: {
      auto & b = ae.as<BinaryAssertionExpr>();
      bool is_impl =
          (b.op == BinaryAssertionOperator::OverlappedImplication
           || b.op == BinaryAssertionOperator::NonOverlappedImplication
           || b.op == BinaryAssertionOperator::Implies);

      if (is_impl) {
        // A `##k` prefix on the antecedent restricts which cycles a
        // match can even start from (the earliest anchor cycle is
        // k); gate the whole implication so it is vacuously true
        // before cycle k instead of dropping the delay, which would
        // otherwise evaluate the consequent (e.g. a `$past` with no
        // real history yet) at cycles no valid match could reach.
        uint32_t lhs_delay = 0;
        const AssertionExpr * lhs_inner = &b.left;
        if (auto lhs_matched = match_const_delay_seq(b.left)) {
          lhs_delay = lhs_matched->first;
          lhs_inner = lhs_matched->second;
        }
        // A plain-Boolean-reducible antecedent (the common case) is
        // handled directly; a multi-element/first-match/nested-clock
        // sequence antecedent (`a ##1 b |-> ...`,
        // `first_match(seq) |-> ...`) falls back to the general
        // bounded sequence matcher.
        Term lhs = assertion_expr_to_bool(*lhs_inner);
        if (!lhs) lhs = match_exists(*lhs_inner);
        if (!lhs) return Term();

        // Compute the consequent at its anchor cycle (offset by any
        // `##k` on the RHS), then delay the antecedent by that
        // offset using a chain of 1-bit latches so the resulting
        // implication is expressed in the current cycle. A *range*
        // delay on a single-element consequent (`##[m:n] Q`) instead
        // becomes an OR over "Q held i cycles ago" for i spanning the
        // window, anchored at the window's latest cycle (delay +=
        // wmax) -- checking at that cycle is exactly when a violation
        // (the whole window has passed with no match) becomes certain.
        uint32_t delay =
            (b.op == BinaryAssertionOperator::NonOverlappedImplication) ? 1 : 0;
        const AssertionExpr * rhs_inner = &b.right;
        Term rhs;
        if (auto matched = match_const_delay_seq(b.right)) {
          delay += matched->first;
          rhs_inner = matched->second;
          rhs = assertion_expr_to_bool(*rhs_inner);
        } else if (b.right.kind == AssertionExprKind::SequenceConcat
                   && b.right.as<SequenceConcatExpr>().elements.size() == 1) {
          auto & elem = b.right.as<SequenceConcatExpr>().elements[0];
          if (!elem.delay.max) {
            throw PonoException(
                "SystemVerilogEncoder: unbounded sequence delay (##[m:$]) "
                "is not supported");
          }
          uint32_t wmin = elem.delay.min;
          uint32_t wmax = *elem.delay.max;
          if (Term inner = assertion_expr_to_bool(*elem.sequence)) {
            delay += wmax;
            rhs = inner;
            for (uint32_t i = 1; i <= wmax - wmin; ++i) {
              rhs = solver_->make_term(Or, rhs, delay_bool(inner, i));
            }
          }
        } else {
          rhs = assertion_expr_to_bool(*rhs_inner);
        }
        if (!rhs) return Term();

        if (delay > 0) {
          // Materialize the antecedent as a 1-bit BV so the
          // history chain has a value to latch.
          Sort bv1 = solver_->make_sort(BV, 1);
          Term one_bv1 = solver_->make_term(1, bv1);
          Term zero_bv1 = solver_->make_term(0, bv1);
          Term lhs_bv = solver_->make_term(Ite, lhs, one_bv1, zero_bv1);
          Term delayed_bv = make_history_chain(lhs_bv, delay);
          lhs = solver_->make_term(Equal, delayed_bv, one_bv1);
        }
        Term result = solver_->make_term(Implies, lhs, rhs);
        if (lhs_delay > 0) {
          result = solver_->make_term(Or, ltl_before_cycle(lhs_delay), result);
        }
        // `disable iff`: exempt this cycle's check if the disable
        // condition held anywhere in the antecedent-to-consequent
        // shift window, not just at the single cycle the check is
        // anchored at.
        if (Term dw = disable_window(delay)) {
          result = solver_->make_term(Or, dw, result);
        }
        return result;
      }

      Term lhs = assertion_expr_to_bool(b.left);
      Term rhs = assertion_expr_to_bool(b.right);
      if (!lhs || !rhs) return Term();
      switch (b.op) {
        case BinaryAssertionOperator::And:
          return solver_->make_term(And, lhs, rhs);
        case BinaryAssertionOperator::Or:
          return solver_->make_term(Or, lhs, rhs);
        case BinaryAssertionOperator::Iff:
          return solver_->make_term(Equal, lhs, rhs);
        default:
          // Intersect / Throughout / Within / Until* / FollowedBy:
          // sequence operators that span multiple cycles in the
          // general case -- out of scope for the current encoder.
          return Term();
      }
    }

    default:
      // Unary temporal operators (`always`, `s_eventually`, ...),
      // FirstMatch, etc.: not supported.  Caller logs the skipped
      // kind.
      return Term();
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

SystemVerilogEncoder::ResolvedOutputAlias
SystemVerilogEncoder::resolve_output_alias(const slang::ast::Symbol * sym) const
{
  bool aliased = false;
  bool has_range = false;
  uint64_t lo = 0;
  uint64_t hi = 0;
  auto alias_it = port_output_aliases_.find(sym);
  while (alias_it != port_output_aliases_.end()) {
    aliased = true;
    const OutputAliasTarget & tgt = alias_it->second;
    if (tgt.has_range) {
      // `sym` (or the range of it we've already narrowed to) is
      // exactly bits [tgt.lo, tgt.hi] of tgt.sym -- remap our
      // accumulated range into tgt.sym's bit numbering.
      if (has_range) {
        lo = tgt.lo + lo;
        hi = tgt.lo + hi;
      } else {
        lo = tgt.lo;
        hi = tgt.hi;
        has_range = true;
      }
    }
    sym = tgt.sym;
    alias_it = port_output_aliases_.find(sym);
  }
  return { sym, aliased, has_range, lo, hi };
}

bool SystemVerilogEncoder::resolve_wire_on_demand(
    const slang::ast::Symbol * sym)
{
  auto it = wire_drivers_.find(sym);
  if (it == wire_drivers_.end()) return false;
  const WireDriver & drv = it->second;
  const void * stmt_key = drv.ca ? static_cast<const void *>(drv.ca)
                                 : static_cast<const void *>(drv.comb);
  if (processed_drivers_.count(stmt_key)) return true;

  if (!resolving_wires_.insert(sym).second) {
    throw PonoException(
        "SystemVerilogEncoder: combinational loop detected involving '"
        + std::string(sym->name) + "'");
  }
  string saved_prefix = prefix_;
  string saved_parent_prefix = parent_prefix_;
  prefix_ = drv.prefix;
  parent_prefix_ = drv.parent_prefix;

  if (drv.ca) {
    process_continuous_assign_once(*drv.ca);
  } else {
    process_always_comb_once(*drv.comb);
  }

  prefix_ = saved_prefix;
  parent_prefix_ = saved_parent_prefix;
  resolving_wires_.erase(sym);
  return true;
}

Term SystemVerilogEncoder::lookup_symbol(const slang::ast::Symbol * sym)
{
  using namespace slang::ast;

  // Procedural for-loop counter: bound to a per-iteration BV
  // constant for the duration of the unrolling.
  auto lvt = loop_var_terms_.find(sym);
  if (lvt != loop_var_terms_.end()) {
    return lvt->second;
  }

  // If `sym` is a child instance's output-port internal, redirect to
  // the parent-side wire so reads resolve to its term.  This may
  // chain through multiple levels of instantiation (e.g. a
  // grandchild's output port connected straight through an
  // intermediate module's own output port, or one element of an
  // instance array wired to a slice of a parent-side bus), so chase
  // it to the root and pick out the resolved bit range, if any.
  auto resolved = resolve_output_alias(sym);
  sym = resolved.sym;
  auto slice = [&](const Term & t) -> Term {
    if (!resolved.has_range) return t;
    return solver_->make_term(Op(Extract, resolved.hi, resolved.lo), t);
  };

  // Wire being defined in the enclosing always_comb block: return
  // the partial accumulated term so that read-modify-write patterns
  // (e.g. `popcount = popcount + din[i];` inside an unrolled for
  // loop) see the previously-written value.
  auto pending_it = pending_comb_updates_.find(sym);
  if (pending_it != pending_comb_updates_.end()) {
    return slice(pending_it->second);
  }

  auto it = symbol_to_term_.find(sym);
  if (it != symbol_to_term_.end()) {
    return slice(it->second);
  }

  // Not resolved yet -- if `sym` is a wire whose driving continuous
  // assign / always_comb block simply hasn't been walked yet (e.g. it
  // appears later in program order than this read), process it now,
  // out of order, and retry.
  if (resolve_wire_on_demand(sym)) {
    auto pit = pending_comb_updates_.find(sym);
    if (pit != pending_comb_updates_.end()) return slice(pit->second);
    auto sit = symbol_to_term_.find(sym);
    if (sit != symbol_to_term_.end()) return slice(sit->second);
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

  // Enum literal (`IDLE`, `REQ`, ...): slang has already evaluated its
  // value at elaboration time, exactly like a parameter.  Declaring an
  // enum-typed *variable* already worked (its type is integral, so
  // type_to_sort() succeeds); only referencing one of the enum's own
  // named values -- ordinary once the enum is declared -- was missing.
  if (sym->kind == SymbolKind::EnumValue) {
    auto & enum_val = sym->as<EnumValueSymbol>();
    const auto & cv = enum_val.getValue();
    if (!cv.isInteger()) {
      throw PonoException("SystemVerilogEncoder: non-integer enum value '"
                          + string(sym->name) + "'");
    }
    auto val = cv.integer();
    uint64_t width = enum_val.getType().getBitWidth();
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

Term SystemVerilogEncoder::wire_seed_term(const slang::ast::Symbol * sym)
{
  auto it = symbol_to_term_.find(sym);
  if (it != symbol_to_term_.end()) return it->second;
  uint64_t width = sym->as<slang::ast::ValueSymbol>().getType().getBitWidth();
  if (width == 0) width = 1;
  Term iv = fts_.make_inputvar(make_name(string(sym->name)),
                               solver_->make_sort(BV, width));
  symbol_to_term_[sym] = iv;
  return iv;
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

Term SystemVerilogEncoder::replace_bits_dynamic(const Term & base,
                                                const Term & slice,
                                                const Term & idx,
                                                uint64_t elem_w)
{
  uint64_t base_w = base->get_sort()->get_width();
  Sort base_sort = solver_->make_sort(BV, base_w);
  Term idx_ext = resize_to(idx, base_w);
  Term shift_amount = idx_ext;
  if (elem_w != 1) {
    Term elem_w_term = solver_->make_term(elem_w, base_sort);
    shift_amount = solver_->make_term(BVMul, idx_ext, elem_w_term);
  }
  Term elem_ones = solver_->make_term(
      BVNot, solver_->make_term(0, solver_->make_sort(BV, elem_w)));
  Term mask =
      solver_->make_term(BVShl, resize_to(elem_ones, base_w), shift_amount);
  Term slice_shifted =
      solver_->make_term(BVShl, resize_to(slice, base_w), shift_amount);
  Term cleared =
      solver_->make_term(BVAnd, base, solver_->make_term(BVNot, mask));
  return solver_->make_term(
      BVOr, cleared, solver_->make_term(BVAnd, slice_shifted, mask));
}

}  // namespace pono
