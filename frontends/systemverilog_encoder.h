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
 ** Parses and elaborates SystemVerilog designs via slang, then converts
 ** the elaborated representation into Pono's FunctionalTransitionSystem
 ** and properties.
 **
 **/

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "core/fts.h"
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
  /** Construct an SystemVerilogEncoder that parses a SystemVerilog file and
   * populates the given FunctionalTransitionSystem.
   *
   *  Supported SystemVerilog constructs (synthesizable subset):
   *    - Module port declarations (input/output/inout)
   *    - always_ff blocks with non-blocking assignments (state elements)
   *    - always_comb blocks and continuous assignments (combinational logic)
   *    - initial blocks (initial state constraints)
   *    - Basic SVA immediate assertions (converted to properties)
   *    - Bitvector types (logic, bit, reg with packed dimensions)
   *
   *  @param filename path to the SystemVerilog source file
   *  @param fts the transition system to populate
   *  @param filelists paths to SystemVerilog list files (".f" files), each
   *         containing one additional source file path per line ('#' and
   *         "//" lines are comments). Relative paths inside a list file
   *         resolve against that list file's directory. All files named by
   *         `filename` and every `filelists` entry are elaborated together
   *         as a single compilation.
   */
  SystemVerilogEncoder(std::string filename,
                       FunctionalTransitionSystem & fts,
                       const std::vector<std::string> & filelists = {});

  ~SystemVerilogEncoder();

  /** @return the vector of safety properties (negated assertions)
   *  found in the design.
   */
  const smt::TermVec & propvec() const { return propvec_; }

  /** @return the per-property generalized-Büchi justice sets
   *  extracted from temporal (LTL) assertions that are not pure
   *  safety properties.  There is one entry per such assertion; each
   *  entry is the set of justice conditions { j_0, ..., j_k } for that
   *  one property, meaning a counterexample to the original assertion
   *  is an infinite lasso along which every j_i holds infinitely
   *  often (i.e. the conjunction of GF j_i).  Feed a single entry's
   *  TermVec straight into LivenessToSafetyTranslator::translate to
   *  reduce that one property to safety.
   *
   *  Each LTL property also installs its own auxiliary tableau state
   *  (next-step "promise" latches, an init flag, and a per-property
   *  activation latch) and transition constraints directly into the
   *  transition system, so the justice sets are only meaningful
   *  against the system the encoder populated.  Because the
   *  per-property activation latch gates that property's time-0
   *  obligation, distinct LTL properties do not interfere: checking
   *  one property's justice set leaves the others' obligations
   *  vacuous.
   */
  const std::vector<smt::TermVec> & ltl_justice() const { return ltl_justice_; }

 private:
  // ---------- Encoding pipeline ----------

  /** Top-level encoding: parse all source files, elaborate, and walk the
   *  design.
   *  @param filename the primary SystemVerilog source file
   *  @param filelists paths to SystemVerilog list files (".f" files) naming
   *         additional source files to parse alongside `filename`
   */
  void encode(const std::string & filename,
              const std::vector<std::string> & filelists);

  /** Process a top-level module instance. */
  void process_module(const slang::ast::InstanceSymbol & inst);

  /** First pass: declare all variables (state vars, inputs, wires).
   *  Walks ports and internal variable declarations.
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

  /** Process an always_comb block to extract combinational definitions. */
  void process_always_comb(const slang::ast::ProceduralBlockSymbol & proc);

  /** Process an initial block to extract initial-state constraints. */
  void process_initial(const slang::ast::ProceduralBlockSymbol & proc);

  /** Process a continuous assignment (assign statement). */
  void process_continuous_assign(const slang::ast::ContinuousAssignSymbol & ca);

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
    NEXT_STATE,     ///< Inside always_ff: build next-state functions
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

  // ---------- Type conversion ----------

  /** Convert a slang type to an SMT sort.
   *  @param type the slang type
   *  @return the corresponding SMT sort (BV or BOOL)
   */
  smt::Sort type_to_sort(const slang::ast::Type & type);

  // ---------- Helpers ----------

  /** Look up the SMT term for a slang symbol.
   *  @param sym pointer to the slang symbol
   *  @return the SMT term, or throws if not found
   */
  smt::Term lookup_symbol(const slang::ast::Symbol * sym);

  /** Result of chasing an output-port alias chain to its root: the
   *  root symbol, whether any hop was taken at all, and -- when the
   *  chain passed through a bus-element connection (e.g. one element
   *  of an instance array wired to a slice of a parent-side array) --
   *  the bit range within the root that this symbol corresponds to.
   */
  struct ResolvedOutputAlias
  {
    const slang::ast::Symbol * sym;
    bool aliased;
    bool has_range;
    uint64_t lo;
    uint64_t hi;
  };

  /** Chase port_output_aliases_ transitively (an output port that is
   *  itself connected to an outer output port, e.g. through two
   *  levels of instantiation, or through a bus-element connection of
   *  an instance array) to the final root symbol and bit range.
   *  @param sym the symbol to resolve
   */
  ResolvedOutputAlias resolve_output_alias(
      const slang::ast::Symbol * sym) const;

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

  /** Ensure a term has the expected bit-width, zero-extending or
   *  truncating as needed.
   *  @param t the term to resize
   *  @param target_width the desired width
   *  @return the resized term
   */
  smt::Term resize_to(const smt::Term & t, uint64_t target_width);

  /** Build a partial-write term: take the full-width `base` and
   *  return a term equal to `base` everywhere except bits
   *  [lo .. hi], which take their values from `slice`.  A full-width
   *  write (lo == 0 && hi == width(base)-1) is just `slice`.
   */
  smt::Term replace_bits(const smt::Term & base,
                         const smt::Term & slice,
                         uint64_t lo,
                         uint64_t hi);

  /** Build a partial-write term like replace_bits(), but for a
   *  runtime-variable element index (`arr[idx] = slice` where `idx`
   *  is not a compile-time constant): shifts an `elem_w`-wide window
   *  of ones into position `idx * elem_w` and uses it to mask the
   *  shifted-into-position `slice` into `base`, leaving every other
   *  element unchanged.  Mirrors the shift+extract technique already
   *  used for dynamic-index *reads* in expr_to_term().
   */
  smt::Term replace_bits_dynamic(const smt::Term & base,
                                 const smt::Term & slice,
                                 const smt::Term & idx,
                                 uint64_t elem_w);

  /** Handle `base[idx] = rhs` (nonblocking or blocking) when `idx` is
   *  not a compile-time constant, so resolve_lvalue() can't produce a
   *  static bit range.  Only a direct select on a plain variable base
   *  is supported (no nested dynamic selects).  A no-op (matching the
   *  prior silent-skip behavior) if the base isn't a plain variable
   *  or has no current term yet.
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
   *  blocks (via pre_scan_always_ff()) and every descendant instance's
   *  body, so every register anywhere in the design tree is known
   *  before any instance's variables are declared -- otherwise a
   *  sibling instance visited earlier in source order (e.g. an
   *  `interface` instance whose members are actually driven by a
   *  later sibling's always_ff through a hierarchical/interface-port
   *  reference) would get its members wrongly declared as free
   *  inputs.
   *  @param body the instance body to scan (and recurse from)
   */
  void pre_scan_state_vars(const slang::ast::InstanceBodySymbol & body);

  /** Pre-scan a combinational always_comb body to identify blocking
   *  assignment targets as combinational wire symbols.  Each wire
   *  target found is recorded (with the current prefix_ /
   *  parent_prefix_) as being driven by `proc`, so a read of it that
   *  occurs before `proc` is naturally walked can force it to be
   *  processed on demand -- see resolve_wire_on_demand().
   *  @param body the statement body of the always_comb block
   *  @param proc the enclosing always_comb block symbol
   */
  void pre_scan_always_comb(const slang::ast::Statement & body,
                            const slang::ast::ProceduralBlockSymbol & proc);

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
   *  @param scope the scope whose members should be visited
   *  @param fn   callback invoked once per concrete member symbol
   */
  template <typename Fn>
  void walk_members(const slang::ast::Scope & scope, Fn && fn);

  /** Lazily construct and return the encoder's slang EvalContext.
   *  Used to evaluate compile-time-constant loop bounds and step
   *  expressions when unrolling procedural `for` loops.
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
   *  (`Simple`, `SequenceConcat` with per-element `[m:n]` delay ranges,
   *  `FirstMatch`, `Clocking` -- unwrapped/ignored per this file's
   *  multiclock design decision), returns a vector indexed by relative
   *  offset `L` where entry `L` is a Term true iff the sequence
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
   *  Throws for a sequence shape this doesn't model (its own leading
   *  repetition, or a Binary intersect/within/throughout as the
   *  outermost sequence) rather than guessing.
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

  /** Build a chain of `n` 1-cycle latch state vars that track
   *  `value` over n cycles, returning a Term equal to the value
   *  from `n` cycles ago.  Each latch is initialised to 0 and its
   *  next-state is the previous link in the chain.  Used by both
   *  `$past(...)` and sequence-delay assertion expressions.
   *  @param value the current-cycle value to delay
   *  @param n     the number of cycles of delay (>= 1)
   *  @return a Term holding `value` from `n` cycles ago
   */
  smt::Term make_history_chain(const smt::Term & value, uint32_t n);

  /** Like make_history_chain(), but for a Bool-sorted (not BV-sorted)
   *  `cond` -- some solvers (e.g. Bitwuzla) can't build a state
   *  variable or a zero constant of sort Bool, so this wraps `cond`
   *  into a BV1 flag before delaying it and unwraps back to Bool
   *  afterward, matching the pattern already used for the `|=>`/
   *  `|-> ##N` antecedent delay.
   *  @param cond the current-cycle Bool-sorted condition to delay
   *  @param n    the number of cycles of delay (>= 1)
   *  @return a Bool-sorted Term equal to `cond` from `n` cycles ago
   */
  smt::Term delay_bool(const smt::Term & cond, uint32_t n);

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
   *  (see `ltl_make_*`) via `assign_next` plus a current-cycle
   *  consistency constraint, and every strong-eventuality operator
   *  (F / strong-until) appends its discharge condition to `justice`.
   *
   *  Returns a null Term when the property uses an operator the
   *  tableau does not model (sequence intersect/throughout/within/
   *  followed-by, etc.); the caller then skips the assertion.
   */
  smt::Term ltl_to_sat(const slang::ast::AssertionExpr & ae,
                       bool neg,
                       smt::TermVec & justice);

  /** Tableau gadget builders.  Each takes the already-compiled
   *  `sat(...)` Boolean term(s) of the operand(s) and returns the
   *  `sat(...)` term of the composite temporal formula, installing the
   *  necessary "promise" latch and consistency constraint into the
   *  transition system.  The F / U builders also push their
   *  eventuality-discharge term onto `justice`.
   */
  smt::Term ltl_make_X(const smt::Term & phi);
  smt::Term ltl_make_G(const smt::Term & phi);
  smt::Term ltl_make_F(const smt::Term & phi, smt::TermVec & justice);
  smt::Term ltl_make_R(const smt::Term & a, const smt::Term & b);
  smt::Term ltl_make_U(const smt::Term & a,
                       const smt::Term & b,
                       smt::TermVec & justice);

  /** Lazily create the shared "first cycle" flag: a 1-bit state var
   *  that is 1 in the initial state and 0 forever after.  Returned as
   *  a Boolean term (`flag == 1`).  Used to gate each LTL property's
   *  time-0 obligation so it is only asserted at the start of the
   *  trace. */
  smt::Term ltl_init_flag();

  /** Returns a Boolean term true iff the current cycle is one of the
   *  first `k` cycles of the trace (cycle index 0..k-1).  Used to
   *  exempt a `##k seq |-> ...` property from being checked before
   *  cycle `k`, the earliest cycle any `##k` match can end at (a
   *  naive per-cycle check would otherwise wrongly evaluate the
   *  consequent -- e.g. a `$past` with no real history yet -- at
   *  cycles no valid match could ever reach).
   *  @param k number of leading cycles to flag (k >= 1)
   */
  smt::Term ltl_before_cycle(uint32_t k);

  /** Returns a Boolean term true iff `current_disable_cond_` holds at
   *  the current cycle or at any of the `window` cycles before it.
   *  Used to extend a `disable iff` condition sampled at an
   *  antecedent's trigger cycle across the cycles a delayed (`##k`)
   *  consequent shifts the check over, so the whole match window is
   *  exempted, not just its last cycle.  Returns a null Term when
   *  `current_disable_cond_` is null (no disable iff in scope).
   *  @param window number of cycles before the current one to OR in
   */
  smt::Term disable_window(uint32_t window);

  // ---------- Data members ----------

  FunctionalTransitionSystem & fts_;
  const smt::SmtSolver & solver_;

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
  struct OutputAliasTarget
  {
    const slang::ast::Symbol * sym;
    bool has_range = false;
    uint64_t lo = 0;
    uint64_t hi = 0;
  };

  // Map from a child instance's port internal symbol to the parent-side
  // variable (and, for a bus-element connection, bit range) that the
  // port connects to.  Populated while processing a child instance and
  // consulted when resolving the LHS of a continuous-assign or
  // always_comb statement so writes to a child output port redirect to
  // the parent-side wire.
  std::unordered_map<const slang::ast::Symbol *, OutputAliasTarget>
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
  // often.  See ltl_justice() and ltl_to_sat().
  std::vector<smt::TermVec> ltl_justice_;

  // Cached Boolean term (`flag == 1`) for the shared LTL "first
  // cycle" state var, created on demand by ltl_init_flag().  Null
  // until the first LTL property is encoded.
  smt::Term ltl_init_flag_;

  // Hierarchical name prefix for the current module.
  std::string prefix_;

  // Hierarchical name prefix for the *parent* of the module currently
  // being processed (i.e., the prefix of the scope where output-port
  // aliases live).  Updated as we recurse into / out of child
  // instances so wires redirected via port_output_aliases_ get named
  // in their owning scope rather than the child's.
  std::string parent_prefix_;

  // Procedural-loop counter bindings: when unrolling a `for` loop,
  // each loop variable's slang VariableSymbol is mapped to a BV
  // term representing its current iteration value.  Consulted by
  // lookup_symbol so the loop body's references resolve to the
  // right constant in each unrolled iteration.
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

  // Monotonic counter used to mint unique state-var names for the
  // hidden latches introduced by `$past`, `|=>`, and `|-> ##N`.
  uint32_t latch_counter_ = 0;

  // Scope of the module instance body currently being processed by
  // process_assignments().  Used to look up that module's `default
  // disable iff` condition (Compilation::getDefaultDisable walks up
  // the scope's parent chain).  Saved/restored around recursion into
  // child instances.
  const slang::ast::Scope * current_scope_ = nullptr;

  // The `disable iff` condition (explicit on the current assert
  // statement, or the enclosing module's `default disable iff`) as a
  // Boolean SMT term, or null if none applies.  Set just before
  // compiling one assertion's property expression and read by
  // assertion_expr_to_bool()/ltl_to_sat() via disable_window(), so it
  // need not be threaded through every recursive call.
  smt::Term current_disable_cond_;
};

}  // namespace pono
