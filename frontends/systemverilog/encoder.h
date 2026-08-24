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
 *                         canonicalization, loop-control-signal type) shared
 *                         across the files below.
 *   - prescan.cpp:        the pre-scan pass that classifies wires vs. state
 *                         vars before declaration.
 *   - declare.cpp:        the variable-declaration pass.
 *   - instance.cpp:       continuous assigns, always_comb/ff/initial blocks,
 *                         and per-instance (child module) processing.
 *   - statement.cpp:      the process_statement() procedural statement encoder.
 *   - expr.cpp:           the expr_to_term() expression encoder.
 *   - sva.cpp:            SVA/LTL AST-walking and dispatch
 * (assertion_expr_to_bool, offsets_ending_now, ltl_to_sat, and friends) --
 *                         needs expr_to_term(), so it stays here rather than
 *                         moving into tableau.{h,cpp}.
 *   - tableau.h/.cpp:     the Tableau class -- the SVA/LTL tableau's pure
 *                         latch-building primitives, with no dependency on
 *                         the rest of this class. Owned by
 *                         SystemVerilogEncoder as tableau_ and called from
 *                         sva.cpp.
 *   - terms.cpp:          low-level bit/term helpers (type-to-sort,
 *                         resize/replace-bits, symbol lookup).
 */

#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "core/fts.h"
#include "frontends/systemverilog/bit_utils.h"
#include "frontends/systemverilog/tableau.h"
#include "smt-switch/smt.h"

// Forward declarations for slang types to avoid exposing slang headers
// to all translation units.
namespace slang::ast {
class Compilation;
class Type;
class Expression;
class Statement;
class Symbol;
class ValueSymbol;
class Scope;
class InstanceSymbol;
class InstanceBodySymbol;
class ProceduralBlockSymbol;
class ContinuousAssignSymbol;
class PortSymbol;
class EvalContext;
class AssertionExpr;
class ElementSelectExpression;
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

  /** First pass: declare state vars and inputs.  Wires are skipped --
   *  they get their term assigned later during combinational-assignment
   *  processing.  Walks ports and internal variable declarations.
   */
  void declare_variables(const slang::ast::InstanceBodySymbol & body);

  /** Declare just the internal (non-port) variables of `body`.  Used
   *  when descending into a child instance, whose ports have already
   *  been bound through the port-connection map.
   */
  void declare_variables_internal(const slang::ast::InstanceBodySymbol & body);

  /** Declare a single port as an input or output variable. */
  void process_port(const slang::ast::PortSymbol & port);

  /** Second pass: process all behavioral and structural assignments.
   *  Walks always blocks, continuous assigns, and initial blocks.
   */
  void process_assignments(const slang::ast::InstanceBodySymbol & body);

  /** Process an always_ff block to extract next-state update functions. */
  void process_always_ff(const slang::ast::ProceduralBlockSymbol & proc);

  /** Shared implementation of process_always_ff(): process `body` with
   *  StmtContext::NEXT_STATE and commit every accumulated
   *  pending_next_updates_ entry via assign_next(). Blocking- vs.
   *  nonblocking-assignment syntax makes no difference to this
   *  encoding -- only the StmtContext does -- so this is also reused
   *  directly for `always_latch` (whose writes use blocking `=`, but
   *  need exactly the same "holds previous value when not written"
   *  register semantics) and for a `forever @(...) ...` loop
   *  recognized as a legacy structural spelling of `always_ff`.
   *  @param body the statement body to process (typically a Timed
   *              statement wrapping the event-controlled block)
   */
  void process_next_state_body(const slang::ast::Statement & body);

  /** Process an always_comb block to extract combinational definitions. */
  void process_always_comb(const slang::ast::ProceduralBlockSymbol & proc);

  /** Process an initial block to extract initial-state constraints. */
  void process_initial(const slang::ast::ProceduralBlockSymbol & proc);

  /** Process a continuous assignment (assign statement). */
  void process_continuous_assign(const slang::ast::ContinuousAssignSymbol & ca);

  /** Process one write from a continuous assignment: `lhs_expr` is
   *  either the assignment's whole LHS (the common case), or one
   *  operand of a concatenation-target LHS (`assign {hi, lo} = ...`),
   *  in which case `rhs` is already sliced down to that operand's own
   *  width.
   *  @param lhs_expr the (non-concatenation) lvalue expression being
   *                  written
   *  @param rhs the value to write, resized to lhs_expr's width
   *  @param rhs_signed whether to sign-extend (rather than zero-
   *         extend) `rhs` if it needs to grow to fit lhs_expr's width
   *         -- the caller passes the *source* expression's own
   *         signedness (false for a concatenation-target slice, which
   *         is positional bit-splicing, not a numeric value)
   */
  void process_continuous_assign_operand(
      const slang::ast::Expression & lhs_expr,
      const smt::Term & rhs,
      bool rhs_signed);

  /** Process an always_comb block, or a continuous assign, unless it
   *  has already been processed (either by the normal walk reaching
   *  it, or because some earlier read forced it via
   *  resolve_wire_on_demand()).  Every call site that would otherwise
   *  call process_always_comb() / process_continuous_assign() directly
   *  goes through these instead, so a driver is never processed twice.
   */
  void process_always_comb_once(const slang::ast::ProceduralBlockSymbol & proc);
  void process_continuous_assign_once(
      const slang::ast::ContinuousAssignSymbol & ca);

  /** If `sym` is a wire whose driving continuous assign / always_comb
   *  block lives in a scope already walked (so its home prefix_ /
   *  parent_prefix_ are known) but hasn't been processed yet -- e.g.
   *  a same-scope wire whose `assign` appears later in program order
   *  than the code that reads it -- process that driver now, out of
   *  order, so lookup_symbol() can retry.  Detects and throws on a
   *  genuine combinational cycle.
   *  @return true if a driver was found (and is now processed, or
   *          already had been); false if `sym` isn't a locally-driven
   *          wire this mechanism knows about.
   */
  bool resolve_wire_on_demand(const slang::ast::Symbol * sym);

  // ---------- Statement processing ----------

  /** Context for statement processing: whether we are building next-state
   *  updates (always_ff), combinational definitions (always_comb), or
   *  initial constraints.
   */
  enum class StmtContext
  {
    NEXT_STATE,     ///< Inside always_ff (also always_latch and a legacy
                    ///< `forever @(...)` spelling of always_ff): build
                    ///< next-state functions
    COMBINATIONAL,  ///< Inside always_comb: build combinational definitions
    INITIAL         ///< Inside initial: build init constraints
  };

  /** Recursively process a statement, extracting assignments.
   *  @param stmt the slang statement to process
   *  @param ctx  what kind of block we are in
   *  @param condition accumulated path condition (for if/case nesting)
   */
  void process_statement(const slang::ast::Statement & stmt,
                         StmtContext ctx,
                         const smt::Term & condition);

  /** Re-derive `loop_var_terms_[&sym]` from `sym`'s current constant
   *  value in eval_ctx() (after a for-loop step, a while/repeat/foreach
   *  iteration, or a plain assignment to a compile-time-unrolled local).
   *  Throws if `sym` isn't a currently-bound integer local. */
  void refresh_loop_var_term(const slang::ast::ValueSymbol & sym);

  // ---------- Expression conversion ----------

  /** Convert a slang expression to an SMT term.
   *  Handles operators, literals, variable references, concatenation,
   *  bit-selects, ternary, etc.
   *  @param expr the slang expression
   *  @return the corresponding SMT term
   */
  smt::Term expr_to_term(const slang::ast::Expression & expr);

  // ---------- Helpers ----------

  /** Look up the SMT term for a slang symbol.
   *  @param sym pointer to the slang symbol
   *  @return the SMT term, or throws if not found
   */
  smt::Term lookup_symbol(const slang::ast::Symbol * sym);

  /** One piece of a (possibly split) output-port alias resolution:
   *  bits [rhs_lo, rhs_hi] of the caller's original write/read window
   *  (0-indexed from the window's own low bit) correspond to bits
   *  [target_lo, target_hi] of `sym`, a final symbol that is itself
   *  guaranteed *not* to be a further output-port alias.
   */
  struct ResolvedAliasPiece
  {
    const slang::ast::Symbol * sym;
    uint64_t target_lo, target_hi;
    uint64_t rhs_lo, rhs_hi;
  };

  /** Chase port_output_aliases_ transitively for bits [lo, hi] of
   *  `sym` (an output port that is itself connected to an outer
   *  output port, e.g. through two levels of instantiation, through a
   *  bus-element connection of an instance array, or through a
   *  concatenation-target connection splitting `sym`'s bits across
   *  several parent-side signals) to one or more final (non-aliased)
   *  pieces. A plain, non-aliased `sym` resolves to a single piece
   *  covering itself. Piece order is unspecified; a caller that needs
   *  to reassemble the whole value should sort by rhs_lo.
   *  @param sym the symbol to resolve
   *  @param lo  the low bit of the window to resolve (relative to sym)
   *  @param hi  the high bit of the window to resolve (relative to sym)
   *  @param rhs_base the rhs-window bit corresponding to `lo` --
   *              callers should omit this (recursive calls set it)
   */
  std::vector<ResolvedAliasPiece> resolve_output_alias_pieces(
      const slang::ast::Symbol * sym,
      uint64_t lo,
      uint64_t hi,
      uint64_t rhs_base = 0) const;

  /** Build the hierarchical name for a symbol.
   *  @param name the local name
   *  @return the fully-qualified name with prefix
   */
  std::string make_name(const std::string & name) const;

  /** Return the existing term for a wire-classified symbol, or --
   *  if this is the first write it has ever received -- create a
   *  fresh free variable sized to its declared bit width to serve as
   *  the splice base for a partial write (e.g. one element of a
   *  multi-driver bus written piecemeal across several instance-array
   *  elements).  The placeholder's own value never surfaces once
   *  every bit of the symbol has been written by some partial write.
   */
  smt::Term wire_seed_term(const slang::ast::Symbol * sym);

  /** Handle `base[idx] = rhs` (nonblocking or blocking) when `idx` is
   *  not a compile-time constant, so resolve_lvalue() can't produce a
   *  static bit range.  Only a direct select on a plain variable base
   *  is supported (no nested dynamic selects).  A no-op if the base
   *  isn't a plain variable, if it resolves through
   *  port_output_aliases_ to anything other than a single whole-symbol
   *  alias, if it has no current term yet, or if `ctx` doesn't apply
   *  to it (a COMBINATIONAL write to a non-wire symbol, or any
   *  INITIAL write, isn't needed by any currently-supported
   *  construct).
   *  @param sel  the dynamic element-select LHS expression
   *  @param rhs_expr the assignment's right-hand side
   *  @param ctx  which kind of block this assignment is in
   *  @param condition accumulated path condition (for if/case nesting)
   */
  void process_dynamic_element_assign(
      const slang::ast::ElementSelectExpression & sel,
      const slang::ast::Expression & rhs_expr,
      StmtContext ctx,
      const smt::Term & condition);

  /** Pre-scan an always_ff body to identify non-blocking assignment
   *  targets as state variable symbols.
   *  @param body the statement body of the always_ff block
   */
  void pre_scan_always_ff(const slang::ast::Statement & body);

  /** Recursively pre-scan an instance body's own `always_ff`/`always`
   *  blocks (via pre_scan_always_ff()), `always_latch` blocks (via
   *  pre_scan_always_latch()), `initial forever @(...)` blocks (the
   *  legacy always_ff spelling, via pre_scan_always_ff()), and every
   *  descendant instance's body, so every register anywhere in the
   *  design tree is known before any instance's variables are
   *  declared -- otherwise a sibling instance visited earlier in
   *  source order (e.g. an `interface` instance whose members are
   *  actually driven by a later sibling's always_ff through a
   *  hierarchical/interface-port reference) would get its members
   *  wrongly declared as free inputs.
   *  @param body the instance body to scan (and recurse from)
   */
  void pre_scan_state_vars(const slang::ast::InstanceBodySymbol & body);

  /** Pre-scan a combinational always_comb body to identify blocking
   *  assignment targets.  A target written only full-width becomes a
   *  combinational wire symbol; a target written (even partly)
   *  through a bit/range select, or written both full- and
   *  partial-width, becomes a state variable instead (mirroring
   *  pre_scan_always_latch()'s always-a-state-var rule for that
   *  case), so the slice can later be constrained with an
   *  add_constraint rather than needing a state term to macro-
   *  substitute into.  Each wire target found is recorded (with the
   *  current prefix_ / parent_prefix_) as being driven by `proc`, so a
   *  read of it that occurs before `proc` is naturally walked can
   *  force it to be processed on demand -- see
   *  resolve_wire_on_demand().
   *  @param body the statement body of the always_comb block
   *  @param proc the enclosing always_comb block symbol
   */
  void pre_scan_always_comb(const slang::ast::Statement & body,
                            const slang::ast::ProceduralBlockSymbol & proc);

  /** Pre-scan an always_latch body to identify every blocking-
   *  assignment target (full- or partial-width alike) as a state
   *  variable. Unlike always_comb's full-vs-partial wire/state-var
   *  split, an always_latch target is *always* a state variable, even
   *  when a single write covers its whole width, since a latch
   *  implicitly holds its previous value along any path that doesn't
   *  reassign it -- exactly the same "defaults to itself" semantics
   *  process_next_state_body()'s NEXT_STATE writes already provide.
   *  @param body the statement body of the always_latch block
   */
  void pre_scan_always_latch(const slang::ast::Statement & body);

  /** Pre-scan a child instance to identify any parent-side variables
   *  that are driven by the child's output ports; those become wires
   *  in the parent's transition system.  Also recurses into nested
   *  instances.
   *  @param inst the child instance to scan
   */
  void pre_scan_instance(const slang::ast::InstanceSymbol & inst);

  /** Process a child module instance: register its port connections,
   *  then walk its body (continuous assigns, always blocks, nested
   *  instances) so the parent's transition system inherits all of
   *  the child's behavior.
   *  @param inst the child instance to process
   */
  void process_instance(const slang::ast::InstanceSymbol & inst);

  /** Invoke `fn` for every concrete member of `scope`, recursing
   *  through generate-block scopes so the caller sees the unrolled
   *  per-iteration members directly.  Uninstantiated generate
   *  blocks (the unselected arm of a generate-if / generate-case)
   *  are skipped.
   *
   *  Takes a type-erased `std::function` (rather than a template
   *  parameter) specifically so its single definition can live in
   *  exactly one systemverilog_*.cpp file while still being callable,
   *  as an ordinary non-template member function, from every other
   *  systemverilog_*.cpp file -- a template's definition would instead
   *  need to be visible (and separately instantiated) in each of those
   *  translation units, which for this method would also force this
   *  header to fully include (rather than just forward-declare) the
   *  slang types its body depends on.
   *  @param scope the scope whose members should be visited
   *  @param fn   callback invoked once per concrete member symbol
   */
  void walk_members(const slang::ast::Scope & scope,
                    const std::function<void(const slang::ast::Symbol &)> & fn);

  /** Lazily construct and return the encoder's slang EvalContext, used
   *  throughout statement/expression processing to evaluate
   *  compile-time-constant expressions -- not just loop bounds and
   *  step expressions when unrolling procedural loops, but also
   *  things like case-statement selectors/patterns, constant
   *  array/bit-select indices, and other constant-foldable operands
   *  (e.g. a `**` exponent or `$past`'s cycle-count argument).
   */
  slang::ast::EvalContext & eval_ctx();

  /** Compile an SVA AssertionExpr into a Boolean SMT term that holds
   *  iff the assertion passes at the current cycle.  Returns a null
   *  Term when the expression uses an unsupported operator
   *  (e.g. liveness, sequence delays inside arbitrary positions,
   *  etc.); the caller can then skip that assertion.
   *  Non-overlapping implication (`|=>`) and the
   *  `|-> ##N` pattern introduce hidden latch state vars so the
   *  "P held N cycles ago" predicate is current-state-only.
   *  @param ae the assertion expression to compile
   *  @return the boolean term, or a null Term when unsupported
   */
  smt::Term assertion_expr_to_bool(const slang::ast::AssertionExpr & ae);

  /** General bounded sequence matching: given a sequence expression
   *  (`Simple`/`SequenceWithMatch` with a consecutive `[m:n]`
   *  repetition, `SequenceConcat` with per-element `[m:n]` delay
   *  ranges, `FirstMatch`, `Clocking` -- unwrapped/ignored per this
   *  file's multiclock design decision -- or a `Binary`
   *  intersect/within/throughout composition of two such sequences),
   *  returns a vector indexed by relative offset `L` where entry `L`
   *  is a Term true iff the sequence
   *  completes a match at the *current* cycle, having started `L`
   *  cycles earlier. A null entry means that offset is structurally
   *  unreachable. Returns an empty vector for sequence shapes this
   *  primitive doesn't (yet) model -- the caller should treat that the
   *  same as an unsupported construct.
   *
   *  Scoped to statically-bounded sequences: an unbounded (`[*]`,
   *  `[+]`, `[*n:$]`) or nonconsecutive/goto repetition, or an
   *  unbounded inter-element delay (`##[m:$]`), throws a clear
   *  PonoException rather than silently mismodeling or dropping it --
   *  this is a permanent architectural boundary of the encoder's
   *  compile-time-bounded model, not a "not implemented yet" gap.
   *  @param seq the sequence expression to match
   *  @return offsets indexed by relative start-to-end span
   */
  smt::TermVec offsets_ending_now(const slang::ast::AssertionExpr & seq);

  /** Convenience wrapper over offsets_ending_now(): ORs together every
   *  reachable offset, i.e. "does `seq` complete a match at the
   *  current cycle, regardless of how long it took". Returns a null
   *  Term if offsets_ending_now() returns no reachable offsets at all
   *  (an unsupported sequence shape) so callers can fall back to their
   *  existing unsupported-construct handling.
   *  @param seq the sequence expression to match
   *  @return the "matches now" Boolean term, or a null Term
   */
  smt::Term match_exists(const slang::ast::AssertionExpr & seq);

  /** The Boolean condition of a sequence's own leading element --
   *  "has an attempt to match `seq` just begun" -- used by
   *  weak_seq_bool() to detect when an in-progress match attempt has
   *  definitely failed. Recurses through FirstMatch/Clocking (like
   *  offsets_ending_now()) and into a SequenceConcat's first element.
   *  Throws for any other sequence shape (its own leading repetition,
   *  a `SequenceWithMatch`, or a `Binary` intersect/within/throughout
   *  as the outermost sequence) rather than guessing.
   *  @param seq the sequence expression to inspect
   *  @return the Boolean "an attempt just started" term
   */
  smt::Term leading_condition(const slang::ast::AssertionExpr & seq);

  /** Builds the `weak(seq)` Boolean safety condition: `seq` carries no
   *  obligation to ever match, but if an attempt began exactly
   *  `S = offsets_ending_now(seq).size() - 1` cycles ago (`S` being
   *  the sequence's own maximum span -- the last possible chance for
   *  that attempt to complete) and no completion happened anywhere in
   *  the intervening window, that attempt has definitely failed.
   *  Checked at every cycle, this covers every possible attempt start
   *  point exactly once, `S` cycles after it began. Returns a null
   *  Term if `seq`'s shape isn't modeled by offsets_ending_now().
   *  @param seq the sequence expression `weak(...)` wraps
   *  @return the Boolean safety term, or a null Term
   */
  smt::Term weak_seq_bool(const slang::ast::AssertionExpr & seq);

  /** General symbolic-tableau translation of an SVA property into the
   *  Boolean SMT term `sat(psi)` that holds at a cycle iff the
   *  (possibly negated) property `psi` holds starting from that cycle,
   *  where `psi` is `ae` when `neg` is false and `!ae` when `neg` is
   *  true.  Negation is pushed through the operators on the fly (so
   *  the gadgets built always correspond to the operators of `psi` in
   *  negation-normal form) -- this keeps the eventuality-fairness
   *  conditions correct regardless of the surrounding polarity.
   *
   *  Each temporal operator instantiates a one-step "promise" latch
   *  (see `tableau_`'s `make_X/G/F/R/U`) via `assign_next`
   *  plus a current-cycle consistency constraint, and every
   *  strong-eventuality operator (F / strong-until) appends its
   *  discharge condition to `justice`.
   *
   *  Returns a null Term when the property uses an operator the
   *  tableau does not model (sequence intersect/throughout/within/
   *  followed-by, etc.); the caller then skips the assertion.
   */
  smt::Term ltl_to_sat(const slang::ast::AssertionExpr & ae,
                       bool neg,
                       smt::TermVec & justice);

  // ---------- Data members ----------

  FunctionalTransitionSystem & fts_;
  const smt::SmtSolver & solver_;

  // The SVA/LTL tableau's pure latch-building primitives (make_X/G/F/R/U,
  // make_history_chain, delay_bool, init_flag, before_cycle,
  // disable_window). Holds no reference back to this class -- called from
  // the AST-walking methods above (assertion_expr_to_bool(), ltl_to_sat(),
  // and friends, in sva.cpp) with explicit Term/count/prefix_ arguments.
  // See tableau.h.
  Tableau tableau_;

  // Unique pointer to slang compilation (hidden from header users)
  std::unique_ptr<slang::ast::Compilation> compilation_;

  // Map from slang symbol pointer to the corresponding SMT term.
  // For state variables, this maps to the current-state term.
  // For wires (continuous-assign / always_comb targets), this maps to
  // the defining SMT expression (macro substitution).
  std::unordered_map<const slang::ast::Symbol *, smt::Term> symbol_to_term_;

  // Set of symbols that are state variables (assigned in always_ff).
  // Populated during the first pass so that the second pass knows
  // which variables are registers vs. wires.
  std::unordered_set<const slang::ast::Symbol *> state_var_symbols_;

  // Set of symbols that are combinational wires (assigned in
  // continuous_assign or always_comb).  These are macro-substituted:
  // their term in symbol_to_term_ is the defining expression itself,
  // not a fresh state or input variable.
  std::unordered_set<const slang::ast::Symbol *> wire_symbols_;

  // A locally-driven wire's driving statement, plus the prefix_ /
  // parent_prefix_ its home scope was walked with -- everything
  // resolve_wire_on_demand() needs to process that statement out of
  // program order.  Exactly one of {ca, comb} is set.
  struct WireDriver
  {
    const slang::ast::ContinuousAssignSymbol * ca = nullptr;
    const slang::ast::ProceduralBlockSymbol * comb = nullptr;
    std::string prefix;
    std::string parent_prefix;
  };
  std::unordered_map<const slang::ast::Symbol *, WireDriver> wire_drivers_;

  // Driving statements (ContinuousAssignSymbol* / ProceduralBlockSymbol*,
  // type-erased) that have already been processed, whether by the
  // normal walk reaching them or via resolve_wire_on_demand() -- so
  // neither path ever processes the same driver twice.
  std::unordered_set<const void *> processed_drivers_;

  // Wires currently being resolved by resolve_wire_on_demand(), to
  // detect a genuine combinational cycle rather than recursing
  // forever.
  std::unordered_set<const slang::ast::Symbol *> resolving_wires_;

  // An output-port alias target: either the whole of `sym` (has_range
  // false), or -- for an instance-array element wired to one slice of
  // a parent-side bus signal -- bits [lo, hi] of `sym`.
  // One segment of an output-port alias: bits [port_lo, port_hi] of
  // the aliased port symbol map affinely onto bits [target_lo,
  // target_hi] of `target`. A plain (non-concatenation) connection
  // has exactly one segment spanning the port's whole width; a
  // concatenation-target connection (`.port({hi, lo})`) has one
  // segment per operand, each covering the corresponding slice of the
  // port's bits.
  struct OutputAliasSegment
  {
    uint64_t port_lo, port_hi;
    const slang::ast::Symbol * target;
    uint64_t target_lo, target_hi;
  };

  // Map from a child instance's port internal symbol to the segment(s)
  // describing which parent-side variable(s) (and bit range(s)) its
  // bits connect to.  Populated while processing a child instance and
  // consulted when resolving the LHS of a continuous-assign or
  // always_comb statement so writes to a child output port redirect to
  // the parent-side wire(s).
  std::unordered_map<const slang::ast::Symbol *,
                     std::vector<OutputAliasSegment>>
      port_output_aliases_;

  // Symbols in pending_comb_updates_ that came from an alias-redirect
  // through port_output_aliases_; their term must be named in
  // parent_prefix_ rather than prefix_ when the always_comb block
  // commits.
  std::unordered_set<const slang::ast::Symbol *> pending_comb_aliased_;

  // For always_ff processing: accumulated conditional next-state updates.
  // Maps state variable term -> conditional next-state expression.
  // After processing the block, these are committed via assign_next().
  std::unordered_map<smt::Term, smt::Term> pending_next_updates_;

  // For always_comb processing: accumulated wire-define expressions
  // keyed by the LHS wire symbol.  After processing the block, the
  // accumulated term is stored in symbol_to_term_.
  std::unordered_map<const slang::ast::Symbol *, smt::Term>
      pending_comb_updates_;

  // Safety properties extracted from SVA assert statements.
  smt::TermVec propvec_;

  // Per-property generalized-Büchi justice sets extracted from
  // temporal (LTL) assertions that are not pure safety.  Each entry
  // is the justice set { j_0, ..., j_k } of one assertion: a
  // counterexample is a lasso along which every j_i holds infinitely
  // often.  See Result::ltl_justice and ltl_to_sat().
  std::vector<smt::TermVec> ltl_justice_;

  // Hierarchical name prefix for the current module.
  std::string prefix_;

  // Hierarchical name prefix for the *parent* of the module currently
  // being processed (i.e., the prefix of the scope where output-port
  // aliases live).  Updated as we recurse into / out of child
  // instances so wires redirected via port_output_aliases_ get named
  // in their owning scope rather than the child's.
  std::string parent_prefix_;

  // Compile-time-constant local bindings: each `for`/`foreach` loop
  // variable is mapped here to a BV term representing its current
  // iteration value while the loop is unrolling, and a plain
  // procedural local declaration (`int x = ...;`) is bound here too,
  // as an immutable per-block constant.  Consulted by lookup_symbol
  // so references resolve to the right constant in each unrolled
  // iteration (or to the local's one bound value).
  std::unordered_map<const slang::ast::Symbol *, smt::Term> loop_var_terms_;

  // Slang evaluation context, lazily constructed on first use to
  // evaluate constant bounds and step expressions during for-loop
  // unrolling.  Owned via unique_ptr because EvalContext is not
  // default-constructible and is only forward-declared here.
  std::unique_ptr<slang::ast::EvalContext> eval_ctx_;

  // Stashed "current value of the LHS" used when expanding compound
  // assignments (`x &= y`, `x += y`, ...).  Slang represents the
  // implicit self-reference in the RHS as an
  // ExpressionKind::LValueReference; expr_to_term returns this term
  // for that case.  Set just before converting a compound RHS and
  // cleared right after, with save/restore for nested contexts.
  smt::Term current_lvalue_term_;

  // Scope of the module instance body currently being processed by
  // process_assignments().  Used to look up that module's `default
  // disable iff` condition (Compilation::getDefaultDisable walks up
  // the scope's parent chain).  Saved/restored around recursion into
  // child instances.
  const slang::ast::Scope * current_scope_ = nullptr;

  // The `disable iff` condition (explicit on the current assert
  // statement, or the enclosing module's `default disable iff`) as a
  // Boolean SMT term, or null if none applies.  Set just before
  // compiling one assertion's property expression, and passed
  // explicitly to tableau_.disable_window() from wherever
  // assertion_expr_to_bool()/ltl_to_sat() needs it, so it need not be
  // threaded through every recursive call in between.
  smt::Term current_disable_cond_;
};

}  // namespace pono
