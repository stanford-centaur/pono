/*!
 * \file symbol_table.h
 * \brief Symbol classification, term binding, and on-demand wire lookup.
 * \author Áron Ricardo Perez-Lopez
 * \date 2026
 * \copyright See the LICENSE file in the top-level source directory.
 *
 * SymbolTable is the SV encoder's symbol table: it classifies every symbol
 * in the design as a state var, a combinational wire, or neither (pass-scan
 * classification, formerly prescan.cpp), binds each declared symbol to its
 * SMT term, and resolves a read of any symbol -- including a wire whose
 * driver hasn't been walked yet, an output-port alias, a loop variable, or
 * a parameter/enum literal -- to a Term (formerly terms.cpp's stateful
 * half).
 *
 * The one genuinely circular dependency in the encoder lives here:
 * resolving a not-yet-processed wire on demand (lookup_symbol() ->
 * resolve_wire_on_demand()) must trigger processing of that wire's driving
 * continuous-assign or always_comb statement, which lives in whatever class
 * walks module/instance bodies. Rather than holding a reference back to
 * that class (which would recreate the same coupling this rearchitecture
 * is trying to remove elsewhere), SymbolTable depends only on the small
 * abstract DriverResolver interface below, implemented by that class and
 * installed via set_driver_resolver().
 *
 * Many of the classification/binding maps are exposed as direct mutable
 * accessors rather than narrow named operations: Declarer, InstanceEncoder,
 * and StatementEncoder each need fine-grained, low-level access to these
 * maps (e.g. binding one symbol's term, checking another's classification)
 * as part of their own per-statement/per-instance logic, so wrapping every
 * such access in its own named method would just relocate the same
 * call-site-level detail behind an extra layer of indirection. This is
 * intentionally a shared, closely-collaborating data structure -- a
 * conventional compiler symbol table -- not an attempt at full information
 * hiding.
 */
#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "core/fts.h"
#include "smt-switch/smt.h"

namespace slang::ast {
class Symbol;
class ContinuousAssignSymbol;
class ProceduralBlockSymbol;
class InstanceSymbol;
class InstanceBodySymbol;
class Scope;
class Statement;
}  // namespace slang::ast

namespace pono {

// A locally-driven wire's driving statement, plus the hierarchical-prefix
// pair its home scope was walked with -- everything resolve_wire_on_demand()
// needs to process that statement out of program order. Exactly one of
// {ca, comb} is set.
struct WireDriver
{
  const slang::ast::ContinuousAssignSymbol * ca = nullptr;
  const slang::ast::ProceduralBlockSymbol * comb = nullptr;
  std::string prefix;
  std::string parent_prefix;
};

// One segment of an output-port alias: bits [port_lo, port_hi] of the
// aliased port symbol map affinely onto bits [target_lo, target_hi] of
// `target`. A plain (non-concatenation) connection has exactly one segment
// spanning the port's whole width; a concatenation-target connection
// (`.port({hi, lo})`) has one segment per operand, each covering the
// corresponding slice of the port's bits.
struct OutputAliasSegment
{
  uint64_t port_lo, port_hi;
  const slang::ast::Symbol * target;
  uint64_t target_lo, target_hi;
};

// One piece of a (possibly split) output-port alias resolution: bits
// [rhs_lo, rhs_hi] of the caller's original write/read window (0-indexed
// from the window's own low bit) correspond to bits [target_lo, target_hi]
// of `sym`, a final symbol that is itself guaranteed *not* to be a further
// output-port alias.
struct ResolvedAliasPiece
{
  const slang::ast::Symbol * sym;
  uint64_t target_lo, target_hi;
  uint64_t rhs_lo, rhs_hi;
};

class SymbolTable
{
 public:
  // See the file-level doc comment for why this exists. `prefix`/
  // `parent_prefix` are the driver's own recorded hierarchical-name
  // scope (captured in WireDriver when the driver was pre-scanned) --
  // the implementer is responsible for switching into that scope
  // (saving/restoring its own current prefix state) around processing
  // the driver, since resolving a wire on demand can happen from a
  // completely different scope than the one that drives it.
  class DriverResolver
  {
   public:
    virtual ~DriverResolver() = default;
    virtual void resolve_continuous_assign(
        const slang::ast::ContinuousAssignSymbol & ca,
        const std::string & prefix,
        const std::string & parent_prefix) = 0;
    virtual void resolve_always_comb(
        const slang::ast::ProceduralBlockSymbol & proc,
        const std::string & prefix,
        const std::string & parent_prefix) = 0;
  };

  SymbolTable(FunctionalTransitionSystem & fts, const smt::SmtSolver & solver);

  void set_driver_resolver(DriverResolver & resolver)
  {
    driver_resolver_ = &resolver;
  }

  // ---------- Classification (formerly prescan.cpp) ----------

  /** Pre-scan: identify state variable symbols by scanning always_ff
   *  blocks for non-blocking assignment targets, before declaring any
   *  variables anywhere in the design -- recurses into every descendant
   *  instance *and checker instance* up front (not just this body's own
   *  direct members) so a sibling instance visited earlier in source
   *  order (e.g. an `interface` instance whose members are actually
   *  driven by a later sibling's always_ff through a hierarchical/
   *  interface-port reference) doesn't get its members wrongly declared
   *  as free inputs before its true driver is discovered. Takes a plain
   *  Scope (not specifically an InstanceBodySymbol) so the same walk
   *  also recurses into a CheckerInstanceBodySymbol, since neither this
   *  function nor walk_members() needs anything InstanceBodySymbol-
   *  specific.
   *  @param body the instance (or checker-instance) body to scan (and
   *         recurse from)
   *  @param prefix the caller's current hierarchical name prefix,
   *         threaded through the recursive member walk
   */
  void pre_scan_state_vars(const slang::ast::Scope & body,
                           std::string & prefix);

  /** Identify every register a clocked block infers; called for every
   *  always_ff/always block in the design.
   *
   *  Non-blocking targets are registers wherever they appear. Blocking
   *  targets are too, but only when `clocked` -- a plain `always` may
   *  be level-sensitive (`always @(*)`), where its blocking writes are
   *  combinational and belong to pre_scan_always_comb() instead. Pass
   *  is_edge_triggered(body) for a plain `always`; always_ff is clocked
   *  by definition.
   */
  void pre_scan_always_ff(const slang::ast::Statement & body,
                          bool clocked,
                          const std::string & prefix);

  /** A fresh unconstrained value standing in for an X, named after
   *  what produced it (`tag` becomes part of the variable's name).
   *  Fresh per occurrence, which is the loosest and so the soundest
   *  reading of "could be anything". */
  /** Symbols an `initial` block writes. One whose value nothing
   *  else drives has to hold it, which can only be decided once
   *  every block has been processed -- see add_initial_only_holds().
   */
  /** Symbols a combinational block assigns on some paths but not
   *  all. Real synthesis infers a latch for one of these, so its
   *  definition becomes a next-state update that falls back to the
   *  symbol's own value rather than a same-cycle equality. */
  std::unordered_set<const slang::ast::Symbol *> & latch_symbols()
  {
    return latch_symbols_;
  }

  std::unordered_set<const slang::ast::Symbol *> & initial_written()
  {
    return initial_written_;
  }

  smt::Term make_unknown_value(const smt::Sort & sort, const std::string & tag);

  /** A read of an unpacked-array cell outside the declared range,
   *  which the LRM gives as X. */
  smt::Term make_out_of_range_value(const smt::Sort & sort);

  /** Give a state variable to each local of `body` that some path
   *  reads before writing (see collect_hold_locals()). Such a local
   *  holds its previous value, which is storage; every other local is
   *  bound to the term its write computes and gets nothing here.
   *
   *  These are declared here rather than by Declarer, whose
   *  walk_members() pass only sees module-scope members -- a local
   *  declared inside a procedural block is a member of that block. */
  void declare_hold_locals(const slang::ast::Statement & body,
                           const std::string & prefix);

  /** Pre-scan an always_latch body to identify every blocking-assignment
   *  target (full- or partial-width alike) as a state variable. Unlike
   *  always_comb's full-vs-partial wire/state-var split, an always_latch
   *  target is *always* a state variable, even when a single write
   *  covers its whole width, since a latch implicitly holds its
   *  previous value along any path that doesn't reassign it.
   */
  void pre_scan_always_latch(const slang::ast::Statement & body);

  /** Pre-scan a combinational always_comb body to identify blocking
   *  assignment targets.  A target written only full-width becomes a
   *  combinational wire symbol; a target written (even partly) through
   *  a bit/range select, or written both full- and partial-width,
   *  becomes a state variable instead (mirroring pre_scan_always_latch()'s
   *  always-a-state-var rule for that case), so the slice can later be
   *  constrained with an add_constraint rather than needing a state term
   *  to macro-substitute into. Each wire target found is recorded (with
   *  the given prefix/parent_prefix) as being driven by `proc`, so a read
   *  of it that occurs before `proc` is naturally walked can force it to
   *  be processed on demand -- see resolve_wire_on_demand() (via
   *  lookup_symbol()).
   *  @param body the statement body of the always_comb block
   *  @param proc the enclosing always_comb block symbol
   *  @param prefix the current hierarchical name prefix
   *  @param parent_prefix the current parent hierarchical name prefix
   */
  void pre_scan_always_comb(const slang::ast::Statement & body,
                            const slang::ast::ProceduralBlockSymbol & proc,
                            const std::string & prefix,
                            const std::string & parent_prefix);

  /** Pre-scan a child instance to identify any parent-side variables
   *  that are driven by the child's output ports; those become wires in
   *  the parent's transition system.  Also recurses into nested
   *  instances.
   *  @param prefix the caller's current hierarchical name prefix,
   *         threaded through the recursive member walk (this function
   *         doesn't build any names itself, but walk_members() needs a
   *         string to use as scratch space -- guaranteed unchanged by
   *         the time this call returns)
   */
  void pre_scan_instance(const slang::ast::InstanceSymbol & inst,
                         std::string & prefix);

  // ---------- Lookup / binding (formerly terms.cpp's stateful half) ----------

  /** Look up the SMT term for a slang symbol, resolving loop-variable
   *  bindings, output-port-alias reconstruction, in-progress
   *  always_comb partial values, on-demand wire resolution (via the
   *  installed DriverResolver), and parameter/enum-literal
   *  materialization, in that order.
   *  @param sym pointer to the slang symbol
   *  @return the SMT term, or throws if not found
   */
  smt::Term lookup_symbol(const slang::ast::Symbol * sym);

  /** Return the existing term for a wire-classified symbol, or -- if
   *  this is the first write it has ever received -- create a fresh
   *  free variable sized to its declared bit width to serve as the
   *  splice base for a partial write.
   */
  smt::Term wire_seed_term(const slang::ast::Symbol * sym,
                           const std::string & prefix);

  /** Chase port_output_aliases_ transitively for bits [lo, hi] of `sym`
   *  to one or more final (non-aliased) pieces. A plain, non-aliased
   *  `sym` resolves to a single piece covering itself. Piece order is
   *  unspecified; a caller that needs to reassemble the whole value
   *  should sort by rhs_lo.
   */
  std::vector<ResolvedAliasPiece> resolve_output_alias_pieces(
      const slang::ast::Symbol * sym,
      uint64_t lo,
      uint64_t hi,
      uint64_t rhs_base = 0) const;

  /** Build the hierarchical name `prefix + "." + name` (or just `name`
   *  if `prefix` is empty).
   */
  std::string make_name(const std::string & prefix,
                        const std::string & name) const;

  // ---------- Direct accessors ----------
  // See the file-level doc comment for why these are exposed directly
  // rather than wrapped in narrower named operations.

  /** Which bit ranges of a shared target have been spliced into it by
   *  an output-port-aliased *register*, one contributor at a time.
   *
   *  A register aliased to only part of its target keeps its own state
   *  var and writes its bits into the target, the way a comb wire
   *  driven from several sibling instances is assembled. Nothing along
   *  that path can tell whether the contributors between them cover
   *  the whole target, so the ranges are recorded here and checked
   *  once every instance has been processed -- an uncovered bit would
   *  otherwise be left free rather than reported.
   */
  std::unordered_map<const slang::ast::Symbol *,
                     std::vector<std::pair<uint64_t, uint64_t>>> &
  spliced_alias_ranges()
  {
    return spliced_alias_ranges_;
  }

  /** The one clock signal and edge this design is allowed to have,
   *  established by whichever property or sampled-value function
   *  names one first, and null until then.
   *
   *  Kept here rather than in AssertionWalker because more than one
   *  place has to agree about it: a property's clocking event and a
   *  `$rose(a, @(posedge clk))`'s belong to the same design, and
   *  checking only the first left the second free to name a clock
   *  the design does not have. `edge` is a slang::ast::EdgeKind held
   *  as an int so this header needs none of slang's enums.
   */
  const slang::ast::Symbol * design_clock_sym() const
  {
    return design_clock_sym_;
  }
  int design_clock_edge() const { return design_clock_edge_; }

  /** Record this clock if none is established yet, and report
   *  whether it agrees with the one that is. */
  bool note_design_clock(const slang::ast::Symbol * sym, int edge)
  {
    if (!design_clock_sym_) {
      design_clock_sym_ = sym;
      design_clock_edge_ = edge;
      return true;
    }
    return sym == design_clock_sym_ && edge == design_clock_edge_;
  }

  std::unordered_map<const slang::ast::Symbol *, smt::Term> & symbol_to_term()
  {
    return symbol_to_term_;
  }
  std::unordered_set<const slang::ast::Symbol *> & state_var_symbols()
  {
    return state_var_symbols_;
  }
  std::unordered_set<const slang::ast::Symbol *> & wire_symbols()
  {
    return wire_symbols_;
  }
  std::unordered_map<const slang::ast::Symbol *, WireDriver> & wire_drivers()
  {
    return wire_drivers_;
  }
  std::unordered_set<const void *> & processed_drivers()
  {
    return processed_drivers_;
  }
  std::unordered_map<const slang::ast::Symbol *,
                     std::vector<OutputAliasSegment>> &
  port_output_aliases()
  {
    return port_output_aliases_;
  }
  std::unordered_set<const slang::ast::Symbol *> & pending_comb_aliased()
  {
    return pending_comb_aliased_;
  }
  std::unordered_map<smt::Term, smt::Term> & pending_next_updates()
  {
    return pending_next_updates_;
  }
  /** Registers written with a *blocking* `=` so far in the clocked
   *  block being walked. A later read of one inside the same block
   *  sees the value just written, not the register's current value --
   *  the only way blocking and non-blocking differ once both are
   *  known to infer a register. Cleared per block alongside
   *  pending_next_updates_. */
  std::unordered_set<const slang::ast::Symbol *> & blocking_next_written()
  {
    return blocking_next_written_;
  }
  std::unordered_map<const slang::ast::Symbol *, smt::Term> &
  pending_comb_updates()
  {
    return pending_comb_updates_;
  }
  std::unordered_map<const slang::ast::Symbol *, smt::Term> & loop_var_terms()
  {
    return loop_var_terms_;
  }

 private:
  /** If `sym` is a wire whose driving continuous assign / always_comb
   *  block lives in a scope already walked (so its home prefix/
   *  parent_prefix are known) but hasn't been processed yet, process
   *  that driver now, out of order, via the installed DriverResolver, so
   *  lookup_symbol() can retry. Detects and throws on a genuine
   *  combinational cycle.
   *  @return true if a driver was found (and is now processed, or
   *          already had been); false if `sym` isn't a locally-driven
   *          wire this mechanism knows about.
   */
  bool resolve_wire_on_demand(const slang::ast::Symbol * sym);

  FunctionalTransitionSystem & fts_;
  const smt::SmtSolver & solver_;
  DriverResolver * driver_resolver_ = nullptr;

  std::unordered_map<const slang::ast::Symbol *, smt::Term> symbol_to_term_;
  std::unordered_set<const slang::ast::Symbol *> state_var_symbols_;
  std::unordered_set<const slang::ast::Symbol *> wire_symbols_;
  std::unordered_map<const slang::ast::Symbol *, WireDriver> wire_drivers_;
  std::unordered_set<const void *> processed_drivers_;
  std::unordered_set<const slang::ast::Symbol *> resolving_wires_;
  std::unordered_map<const slang::ast::Symbol *,
                     std::vector<OutputAliasSegment>>
      port_output_aliases_;
  std::unordered_set<const slang::ast::Symbol *> pending_comb_aliased_;
  std::unordered_map<const slang::ast::Symbol *,
                     std::vector<std::pair<uint64_t, uint64_t>>>
      spliced_alias_ranges_;
  const slang::ast::Symbol * design_clock_sym_ = nullptr;
  int design_clock_edge_ = 0;
  std::unordered_map<smt::Term, smt::Term> pending_next_updates_;
  std::unordered_set<const slang::ast::Symbol *> blocking_next_written_;
  std::unordered_set<const slang::ast::Symbol *> latch_symbols_;

  std::unordered_set<const slang::ast::Symbol *> initial_written_;

  uint64_t unknown_counter_ = 0;
  std::unordered_map<const slang::ast::Symbol *, smt::Term>
      pending_comb_updates_;
  std::unordered_map<const slang::ast::Symbol *, smt::Term> loop_var_terms_;
};

}  // namespace pono
