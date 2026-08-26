/*!
 * \file encoder.h
 * \brief SystemVerilog frontend encoder using the slang library.
 * \author Áron Ricardo Perez-Lopez
 * \date 2026
 * \copyright See the LICENSE file in the top-level source directory.
 *
 * Parses and elaborates SystemVerilog designs via slang, then converts
 * the elaborated representation into Pono's FunctionalTransitionSystem
 * and properties.
 *
 * This class's declaration lives in one file, but its implementation is
 * split by concern across several other files in this directory
 * (frontends/systemverilog/):
 *
 *   - encoder.cpp:        the public encode() factory, construction, slang
 *                         compilation, and the top-level per-module
 *                         encoding dispatch (run(), process_module()).
 *   - ast_helpers.h/.cpp: free-function AST helpers (lvalue resolution, modport
 *                         canonicalization, loop-control-signal type,
 *                         walk_members()) shared across the files below.
 *   - bit_utils.h/.cpp:   free-function pure bit/type helpers (type-to-sort,
 *                         resize/replace-bits), with no dependency on this
 *                         class at all.
 *   - symbol_table.h/.cpp: the SymbolTable class -- classifies wires vs.
 *                         state vars before declaration (formerly
 *                         prescan.cpp) and resolves a symbol to its SMT
 *                         term, including on-demand wire resolution
 *                         (formerly terms.cpp's stateful half). Owned by
 *                         SystemVerilogEncoder as symbol_table_.
 *                         InstanceEncoder implements its
 *                         SymbolTable::DriverResolver interface (see
 *                         symbol_table.h) -- the one genuinely circular
 *                         dependency in this whole rearchitecture.
 *   - expr_encoder.h/.cpp: the ExprEncoder class -- expr_to_term() and its
 *                         shared EvalContext, depending only on SymbolTable
 *                         (to resolve a read) and Tableau (for the sampled-
 *                         value system functions). Owned by
 *                         SystemVerilogEncoder as expr_encoder_.
 *   - declarer.h/.cpp:    the Declarer class -- the variable-declaration
 *                         pass, depending only on SymbolTable. Owned by
 *                         SystemVerilogEncoder as declarer_.
 *   - instance_encoder.h/.cpp: the InstanceEncoder class -- continuous
 *                         assigns, always_comb/ff/initial blocks, and
 *                         per-instance (child module) processing,
 *                         depending on SymbolTable, Declarer,
 *                         StatementEncoder, and ExprEncoder. Owned by
 *                         SystemVerilogEncoder as instance_encoder_.
 *   - statement_encoder.h/.cpp: the StatementEncoder class -- the
 *                         process_statement() procedural statement encoder,
 *                         depending on SymbolTable, ExprEncoder, and
 *                         AssertionWalker. Owned by SystemVerilogEncoder as
 *                         statement_encoder_.
 *   - assertion_walker.h/.cpp: the AssertionWalker class -- concurrent/
 *                         immediate assertion-statement dispatch and SVA/LTL
 *                         AST-walking (assertion_expr_to_bool,
 *                         offsets_ending_now, ltl_to_sat, and friends),
 *                         depending only on ExprEncoder and Tableau. Owned by
 *                         SystemVerilogEncoder as assertion_walker_, which
 *                         also owns the extracted propvec()/ltl_justice().
 *   - tableau.h/.cpp:     the Tableau class -- the SVA/LTL tableau's pure
 *                         latch-building primitives, with no dependency on
 *                         the rest of this class. Owned by
 *                         SystemVerilogEncoder as tableau_ and called from
 *                         assertion_walker.cpp.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "core/fts.h"
#include "frontends/systemverilog/assertion_walker.h"
#include "frontends/systemverilog/declarer.h"
#include "frontends/systemverilog/expr_encoder.h"
#include "frontends/systemverilog/instance_encoder.h"
#include "frontends/systemverilog/statement_encoder.h"
#include "frontends/systemverilog/symbol_table.h"
#include "frontends/systemverilog/tableau.h"
#include "smt-switch/smt.h"

// Forward declarations for slang types to avoid exposing slang headers
// to all translation units.
namespace slang::ast {
class Compilation;
class InstanceSymbol;
}  // namespace slang::ast

namespace pono {

class SystemVerilogEncoder
{
 public:
  /** The properties extracted from a SystemVerilog design by encode(). */
  struct Result
  {
    /** The vector of safety properties (negated assertions) found in
     *  the design.
     */
    smt::TermVec propvec;

    /** The per-property generalized-Büchi justice sets extracted from
     *  temporal (LTL) assertions that are not pure safety properties.
     *  There is one entry per such assertion; each entry is the set
     *  of justice conditions { j_0, ..., j_k } for that one property,
     *  meaning a counterexample to the original assertion is an
     *  infinite lasso along which every j_i holds infinitely often
     *  (i.e. the conjunction of GF j_i).  Feed a single entry's
     *  TermVec straight into LivenessToSafetyTranslator::translate to
     *  reduce that one property to safety.
     *
     *  Each LTL property also installs its own auxiliary tableau
     *  state (next-step "promise" latches and a per-property
     *  activation latch, gated by the one shared "first cycle" init
     *  flag common to every LTL property) and transition constraints
     *  directly into the transition system, so the justice sets are
     *  only meaningful against the system encode() populated.
     *  Because the per-property activation latch gates that
     *  property's time-0 obligation, distinct LTL properties do not
     *  interfere: checking one property's justice set leaves the
     *  others' obligations vacuous.
     */
    std::vector<smt::TermVec> ltl_justice;
  };

  /** Parse a SystemVerilog file, elaborate it, and encode it into the
   *  given FunctionalTransitionSystem.
   *
   *  Supported SystemVerilog constructs (synthesizable subset):
   *    - Module port declarations (input/output/inout)
   *    - always_ff blocks with non-blocking assignments (state elements)
   *    - always_comb blocks and continuous assignments (combinational logic)
   *    - initial blocks (initial state constraints)
   *    - Basic SVA immediate assertions (converted to properties)
   *    - Bitvector types (logic, bit, reg with packed dimensions)
   *
   *  Error handling contract for anything outside that subset -- every
   *  construct this encoder doesn't support falls into exactly one of
   *  two buckets, never a third "silently produces a wrong encoding"
   *  bucket:
   *    - Simulation-only constructs (no possible functional-logic
   *      meaning in a per-cycle model at all: `program`/`checker`
   *      instances, `specify` blocks, `fork`/`join`, `wait`,
   *      `force`/`release`, `defparam`, `bind`, ...) are dropped and
   *      logged as a warning (`logger.log(1, ...)`, visible at
   *      verbosity >= 1) rather than encoded or rejected -- the rest
   *      of the design is still encoded normally.
   *    - Everything else this encoder doesn't (yet) support -- a
   *      mainstream RTL feature not yet implemented, or a malformed/
   *      out-of-bounds use of a feature it does implement -- throws a
   *      PonoException instead of silently dropping or mis-encoding
   *      it. This is enforced structurally, not just by convention:
   *      e.g. resolve_lvalue()'s default case throws for any lvalue
   *      expression shape it has no case for, rather than returning
   *      nullopt and letting the caller silently no-op the write.
   *  See tests/encoders/test_systemverilog_unsupported.cpp for the
   *  ledger of constructs checked against this contract, each verified
   *  empirically (not assumed) to land in the bucket its test expects.
   *
   *  @param filename path to the SystemVerilog source file
   *  @param fts the transition system to populate
   *  @param filelists paths to SystemVerilog list files (".f" files), each
   *         containing one additional source file path per line ('#' and
   *         "//" lines are comments). Relative paths inside a list file
   *         resolve against that list file's directory. All files named by
   *         `filename` and every `filelists` entry are elaborated together
   *         as a single compilation.
   *  @return the safety properties and LTL justice sets found in the design
   */
  static Result encode(std::string filename,
                       FunctionalTransitionSystem & fts,
                       const std::vector<std::string> & filelists = {});

  ~SystemVerilogEncoder();

 private:
  /** Construct a SystemVerilogEncoder bound to `fts` and its solver, doing
   *  no parsing or encoding -- see encode(), the only public entry point,
   *  which is also the only caller of this constructor.
   */
  explicit SystemVerilogEncoder(FunctionalTransitionSystem & fts);

  // ---------- Encoding pipeline ----------

  /** Top-level encoding: parse all source files, elaborate, and walk the
   *  design.  Called exactly once, by encode(), right after construction.
   *  @param filename the primary SystemVerilog source file
   *  @param filelists paths to SystemVerilog list files (".f" files) naming
   *         additional source files to parse alongside `filename`
   */
  void run(const std::string & filename,
           const std::vector<std::string> & filelists);

  /** Process a top-level module instance. */
  void process_module(const slang::ast::InstanceSymbol & inst);

  // ---------- Data members ----------

  FunctionalTransitionSystem & fts_;
  const smt::SmtSolver & solver_;

  // The SVA/LTL tableau's pure latch-building primitives (make_X/G/F/R/U,
  // make_history_chain, delay_bool, init_flag, before_cycle,
  // disable_window). Holds no reference back to this class -- called from
  // the AST-walking methods above (assertion_expr_to_bool(), ltl_to_sat(),
  // and friends, in assertion_walker.cpp) with explicit Term/count/prefix
  // arguments. See tableau.h.
  Tableau tableau_;

  // Unique pointer to slang compilation (hidden from header users)
  std::unique_ptr<slang::ast::Compilation> compilation_;

  // Symbol classification, term binding, and on-demand wire lookup --
  // see symbol_table.h. Holds no reference back to this class;
  // instance_encoder_ implements symbol_table_'s DriverResolver
  // interface (see instance_encoder.h).
  SymbolTable symbol_table_;

  // expr_to_term() and its shared EvalContext -- see expr_encoder.h.
  // Holds no reference back to this class.
  ExprEncoder expr_encoder_;

  // The variable-declaration pass -- see declarer.h. Holds no reference
  // back to this class.
  Declarer declarer_;

  // Assertion-statement dispatch and SVA/LTL AST-walking -- see
  // assertion_walker.h. Holds no reference back to this class; owns the
  // extracted propvec()/ltl_justice() that Result is built from.
  AssertionWalker assertion_walker_;

  // The process_statement() procedural statement encoder -- see
  // statement_encoder.h. Holds no reference back to this class.
  StatementEncoder statement_encoder_;

  // Continuous assigns, always_comb/ff/initial blocks, and per-instance
  // (child module) processing -- see instance_encoder.h. Holds no
  // reference back to this class; implements symbol_table_'s
  // DriverResolver interface for real.
  InstanceEncoder instance_encoder_;
};

}  // namespace pono
