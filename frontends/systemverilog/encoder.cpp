/*!
 * \file encoder.cpp
 * \brief Top-level pipeline: source loading, slang compilation, and dispatch.
 * \author Áron Ricardo Perez-Lopez
 * \date 2026
 * \copyright See the LICENSE file in the top-level source directory.
 *
 * Besides the public encode() factory, run() (the pipeline it delegates to
 * right after constructing the encoder), and process_module(), this file's
 * anonymous namespace parses dot-f list files (rejecting unsupported +/-
 * tool directives) and scans the raw syntax tree for `bind` directives,
 * which have no elaborated Symbol to catch during the normal member walk
 * (walk_members() itself is a free function in ast_helpers.h/.cpp).
 * process_module() runs four ordered passes -- state-variable pre-scan,
 * combinational-wire pre-scan, variable declaration, then assignment
 * processing -- since later passes rely on symbol classifications
 * established earlier.
 */

#include "frontends/systemverilog/encoder.h"

#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "frontends/systemverilog/ast_helpers.h"
#include "slang/ast/ASTContext.h"
#include "slang/ast/Compilation.h"
#include "slang/ast/EvalContext.h"
#include "slang/ast/Expression.h"
#include "slang/ast/Scope.h"
#include "slang/ast/SemanticFacts.h"
#include "slang/ast/expressions/AssignmentExpressions.h"
#include "slang/ast/expressions/MiscExpressions.h"
#include "slang/ast/symbols/BlockSymbols.h"
#include "slang/ast/symbols/CompilationUnitSymbols.h"
#include "slang/ast/symbols/InstanceSymbols.h"
#include "slang/ast/symbols/MemberSymbols.h"
#include "slang/diagnostics/DiagnosticEngine.h"
#include "slang/syntax/AllSyntax.h"
#include "slang/syntax/SyntaxTree.h"
#include "utils/exceptions.h"
#include "utils/logger.h"

using namespace std;

namespace pono {

// ============================================================================
// Construction / destruction
// ============================================================================

SystemVerilogEncoder::SystemVerilogEncoder(FunctionalTransitionSystem & fts)
    : fts_(fts),
      solver_(fts.solver()),
      tableau_(fts_, solver_),
      symbol_table_(fts_, solver_),
      expr_encoder_(symbol_table_, tableau_, solver_),
      declarer_(symbol_table_, fts_, solver_),
      assertion_walker_(expr_encoder_, tableau_, solver_, fts_),
      statement_encoder_(
          symbol_table_, expr_encoder_, assertion_walker_, fts_, solver_),
      instance_encoder_(symbol_table_,
                        declarer_,
                        statement_encoder_,
                        expr_encoder_,
                        fts_,
                        solver_)
{
}

SystemVerilogEncoder::~SystemVerilogEncoder() = default;

SystemVerilogEncoder::Result SystemVerilogEncoder::encode(
    string filename,
    FunctionalTransitionSystem & fts,
    const std::vector<std::string> & filelists)
{
  SystemVerilogEncoder enc(fts);
  enc.run(filename, filelists);
  return Result{ std::move(enc.assertion_walker_.propvec()),
                 std::move(enc.assertion_walker_.ltl_justice()) };
}

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

// Recursively scans `node` for `bind` directive syntax, logging a
// warning for each one found. A `bind` directive injects an
// additional instance into another module's scope from outside its
// source and has no functional-logic representation in this encoder's
// model, but (unlike a program/checker instance or a specify block)
// it has no corresponding elaborated Symbol to catch during the usual
// member walk -- it has to be found in the raw syntax tree instead.
void warn_on_bind_directives(const slang::syntax::SyntaxNode & node)
{
  if (node.kind == slang::syntax::SyntaxKind::BindDirective) {
    logger.log(
        1,
        "SystemVerilogEncoder: ignoring `bind` directive (simulation-only "
        "construct)");
  }
  for (size_t i = 0, n = node.getChildCount(); i < n; i++) {
    if (auto * child = node.childNode(i)) {
      warn_on_bind_directives(*child);
    }
  }
}

}  // namespace

// ============================================================================
// Top-level encoding pipeline
// ============================================================================

void SystemVerilogEncoder::run(const string & filename,
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
  expr_encoder_.bind_compilation(*compilation_);
  instance_encoder_.bind_compilation(*compilation_);
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

  // `bind` directives have no corresponding elaborated Symbol to catch
  // during the usual member walk, so scan the raw syntax for them
  // directly and warn -- see warn_on_bind_directives().
  for (auto & tree : compilation_->getSyntaxTrees()) {
    warn_on_bind_directives(tree->root());
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
  string prefix(inst.name);
  logger.log(1, "SystemVerilogEncoder: processing module {}", prefix);
  // The top-level module has no enclosing scope, so its own
  // parent_prefix (where its own output-port aliases would live) is
  // empty.
  string parent_prefix;

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
  symbol_table_.pre_scan_state_vars(body, prefix);

  // Second pre-pass: identify combinational wire symbols from
  // always_comb blocks, legacy `always` blocks without non-blocking
  // assignments, and continuous-assign LHS values.  Wires are not
  // independent variables -- they will be macro-substituted with their
  // defining expressions, so we must skip declaring them as input vars.
  walk_members(body, prefix, [&](const slang::ast::Symbol & member) {
    if (member.kind == slang::ast::SymbolKind::ProceduralBlock) {
      auto & proc = member.as<slang::ast::ProceduralBlockSymbol>();
      if (proc.procedureKind == slang::ast::ProceduralBlockKind::AlwaysComb) {
        symbol_table_.pre_scan_always_comb(
            proc.getBody(), proc, prefix, parent_prefix);
      } else if (proc.procedureKind
                 == slang::ast::ProceduralBlockKind::Always) {
        // Legacy always: combinational iff it has no non-blocking
        // assignments to identify it as sequential.
        std::unordered_set<const slang::ast::Symbol *> nb_targets;
        collect_nonblocking_targets(proc.getBody(), nb_targets);
        if (nb_targets.empty()) {
          symbol_table_.pre_scan_always_comb(
              proc.getBody(), proc, prefix, parent_prefix);
        }
      }
    } else if (member.kind == slang::ast::SymbolKind::ContinuousAssign) {
      auto & ca = member.as<slang::ast::ContinuousAssignSymbol>();
      auto & ae = ca.getAssignment();
      if (ae.kind == slang::ast::ExpressionKind::Assignment) {
        auto & lhs = ae.as<slang::ast::AssignmentExpression>().left();
        if (lhs.kind == slang::ast::ExpressionKind::NamedValue) {
          auto * sym = &lhs.as<slang::ast::NamedValueExpression>().symbol;
          if (!symbol_table_.state_var_symbols().count(sym)) {
            symbol_table_.wire_symbols().insert(sym);
            symbol_table_.wire_drivers()[sym] = {
              &ca, nullptr, prefix, parent_prefix
            };
          }
        } else if (auto * base = find_lhs_base(lhs)) {
          // Partial-LHS continuous assign (`assign arr[i] = ...`):
          // the base needs to be a state var so process_continuous_assign
          // can constrain the slice via add_constraint.
          if (!symbol_table_.wire_symbols().count(base)) {
            symbol_table_.state_var_symbols().insert(base);
          }
        }
      }
    } else if (member.kind == slang::ast::SymbolKind::Instance) {
      symbol_table_.pre_scan_instance(member.as<slang::ast::InstanceSymbol>(),
                                      prefix);
    }
  });

  // Third pass: declare state vars, inputs, and output-port aliases.
  // Wire symbols are skipped here -- they get their defining term
  // assigned during combinational-assignment processing.
  declarer_.declare_variables(body, prefix);

  // Fourth pass: process behavioral code and continuous assignments.
  instance_encoder_.process_assignments(body, prefix, parent_prefix);
}

}  // namespace pono
