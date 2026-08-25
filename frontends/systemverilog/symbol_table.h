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
   *  instance up front (not just this body's own direct members) so a
   *  sibling instance visited earlier in source order (e.g. an
   *  `interface` instance whose members are actually driven by a later
   *  sibling's always_ff through a hierarchical/interface-port
   *  reference) doesn't get its members wrongly declared as free inputs
   *  before its true driver is discovered.
   *  @param body the instance body to scan (and recurse from)
   *  @param prefix the caller's current hierarchical name prefix,
   *         threaded through the recursive member walk
   */
  void pre_scan_state_vars(const slang::ast::InstanceBodySymbol & body,
                           std::string & prefix);

  /** Thin wrapper around collect_nonblocking_targets(); called for every
   *  always_ff/always block in the design.
   */
  void pre_scan_always_ff(const slang::ast::Statement & body);

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
  std::unordered_map<smt::Term, smt::Term> pending_next_updates_;
  std::unordered_map<const slang::ast::Symbol *, smt::Term>
      pending_comb_updates_;
  std::unordered_map<const slang::ast::Symbol *, smt::Term> loop_var_terms_;
};

}  // namespace pono
