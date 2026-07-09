#ifdef WITH_SLANG

#include <string>
#include <vector>

#include "core/fts.h"
#include "engines/bmc.h"
#include "frontends/systemverilog_encoder.h"
#include "gtest/gtest.h"
#include "modifiers/control_signals.h"
#include "modifiers/liveness_to_safety_translator.h"
#include "modifiers/mod_ts_prop.h"
#include "smt-switch/utils.h"
#include "smt/available_solvers.h"
#include "test_encoder_inputs.h"
#include "utils/exceptions.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class SVUnitTests : public ::testing::Test,
                    public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  static string sv_path(const string & name)
  {
    return string(STRFY(PONO_SRC_DIR)) + "/tests/encoders/inputs/systemverilog/"
           + name;
  }

  // Locate the `rst` input port in the encoded transition system,
  // returning a null Term if the design has no such input.  The
  // encoder names ports as "<top_module>.<port>", so we scan
  // named_terms for any entry ending in ".rst".
  static Term find_reset(const TransitionSystem & ts)
  {
    const string suffix = ".rst";
    for (const auto & [name, term] : ts.named_terms()) {
      if (name.size() >= suffix.size()
          && name.compare(name.size() - suffix.size(), suffix.size(), suffix)
                 == 0) {
        return term;
      }
    }
    return Term();
  }

  // Encode `file` (which must contain a single assert property), run
  // BMC up to `bound` cycles, and check that the result matches
  // `expected`.  When `expected` is FALSE, the resulting witness
  // length must equal `bound` -- i.e., the property fails at exactly
  // cycle `bound` and not earlier.
  //
  // If the design has a top-level `rst` input, the same reset
  // preprocessing that pono.cpp does for --reset is applied: rst is
  // held high for one cycle and the property is guarded with the
  // resulting reset_done term, so each test can pin a deterministic
  // initial state via rst rather than relying on free initial
  // register values.  Designs without an `rst` input (e.g. ones that
  // already pin their state via `initial`, or that are purely
  // combinational) skip this step.
  void check_bmc(const string & file,
                 size_t bound,
                 ProverResult expected = ProverResult::FALSE,
                 const vector<string> & filelists = {})
  {
    SmtSolver s = create_solver(GetParam());
    s->set_opt("incremental", "true");
    s->set_opt("produce-models", "true");
    FunctionalTransitionSystem fts(s);
    SystemVerilogEncoder enc(sv_path(file), fts, filelists);
    ASSERT_EQ(enc.propvec().size(), 1u);
    Term prop_term = enc.propvec()[0];

    TransitionSystem ts = fts;
    if (Term rst = find_reset(ts)) {
      Term reset_done = add_reset_seq(ts, rst, /*reset_bnd=*/1);
      prop_term = ts.solver()->make_term(Implies, reset_done, prop_term);
    }

    // Promote any input vars referenced by the property to state vars,
    // matching the rest of pono.cpp's preprocessing.  SafetyProver
    // only accepts predicates over current-state variables.
    if (!ts.only_curr(prop_term) && ts.no_next(prop_term)) {
      UnorderedTermSet ivs_in_prop;
      get_free_symbolic_consts(prop_term, ivs_in_prop);
      ts = promote_inputvars(ts, ivs_in_prop);
    }

    SafetyProperty prop(ts.solver(), prop_term);
    Bmc bmc(prop, ts, s);
    EXPECT_EQ(bmc.check_until(bound), expected);
    if (expected == ProverResult::FALSE) {
      EXPECT_EQ(bmc.witness_length(), bound);
    }
  }

  // Encode `file` (which must contain a single temporal/LTL
  // assertion), reduce its generalized-Büchi justice set to a safety
  // property via LivenessToSafetyTranslator, then run BMC up to
  // `bound`.  Reset preprocessing is applied if the design has an
  // `rst` input.  The L2S translator adds a "save" input and a few
  // bookkeeping state vars that BMC chooses freely; a FALSE verdict
  // means BMC found a fair lasso violating the LTL property.
  void check_liveness_bmc(const string & file,
                          size_t bound,
                          ProverResult expected = ProverResult::FALSE)
  {
    SmtSolver s = create_solver(GetParam());
    s->set_opt("incremental", "true");
    s->set_opt("produce-models", "true");
    FunctionalTransitionSystem fts(s);
    SystemVerilogEncoder enc(sv_path(file), fts);
    ASSERT_EQ(enc.propvec().size(), 0u);
    ASSERT_EQ(enc.ltl_justice().size(), 1u);
    TermVec justice = enc.ltl_justice()[0];

    TransitionSystem ts = fts;
    if (Term rst = find_reset(ts)) {
      Term reset_done = add_reset_seq(ts, rst, /*reset_bnd=*/1);
      // Pin the accepting lasso to the post-reset region by requiring
      // every justice condition to fire while reset has released.
      for (Term & j : justice) {
        j = ts.solver()->make_term(And, reset_done, j);
      }
    }

    Term safety_term = LivenessToSafetyTranslator{}.translate(ts, justice);

    if (!ts.only_curr(safety_term) && ts.no_next(safety_term)) {
      UnorderedTermSet ivs_in_prop;
      get_free_symbolic_consts(safety_term, ivs_in_prop);
      ts = promote_inputvars(ts, ivs_in_prop);
    }

    SafetyProperty prop(ts.solver(), safety_term);
    Bmc bmc(prop, ts, s);
    EXPECT_EQ(bmc.check_until(bound), expected);
  }
};

// ---------------------------------------------------------------------------
// Block-level constructs
// ---------------------------------------------------------------------------

// Module ports + always_ff with non-blocking assignment + binary +.
// After the one-cycle reset, count starts at 0 and reaches 4 five
// cycles later.
TEST_P(SVUnitTests, EncodeCounter) { check_bmc("counter.sv", 5); }

// `initial` block sets register's initial state.  No rst signal --
// the initial constraint already pins the state, so the property
// fails at cycle 0.
TEST_P(SVUnitTests, InitialBlock) { check_bmc("initial_block.sv", 0); }

// Sized integer literals (hex, binary, signed negative).  After
// reset, sel=1 drives all three registers to the violating literals
// one cycle later.
TEST_P(SVUnitTests, SizedLiterals) { check_bmc("sized_literals.sv", 2); }

// Multiple `assert property` statements yield multiple props.
TEST_P(SVUnitTests, MultipleAssertions)
{
  SmtSolver s = create_solver(GetParam());
  FunctionalTransitionSystem fts(s);
  SystemVerilogEncoder enc(sv_path("multi_assert.sv"), fts);
  EXPECT_EQ(enc.propvec().size(), 3u);
}

// Hierarchical SystemVerilog: a parent module instantiates a child
// that drives one of the parent's wires.  Sub-modules are now
// encoded, so the child's `assign` constrains the wire and BMC
// cannot falsify the assertion at any finite bound.
TEST_P(SVUnitTests, HierarchicalModules)
{
  check_bmc("hierarchical.sv", 5, ProverResult::UNKNOWN);
}

// Hierarchical value reference: the parent reads a child's internal
// register via a dotted name (`c.cnt`) rather than through a port.
// Exercises ExpressionKind::HierarchicalValue resolution.  c.cnt
// hits 5 at cycle 6 (one reset cycle + five increments from 0).
TEST_P(SVUnitTests, HierarchicalValue)
{
  check_bmc("hierarchical_value.sv", 6);
}

// --sv-filelist support: `top` (filelist_top.sv) instantiates `child`,
// which is named only by filelist.f rather than appearing in
// filelist_top.sv itself.  Mirrors HierarchicalModules: if the
// list-file's extra source is elaborated, BMC cannot falsify the
// assertion at any finite bound.
TEST_P(SVUnitTests, Filelist)
{
  check_bmc(
      "filelist_top.sv", 5, ProverResult::UNKNOWN, { sv_path("filelist.f") });
}

// A list file with an unsupported directive line (e.g. "+incdir+...")
// must be rejected rather than silently misread as a filename.
TEST_P(SVUnitTests, FilelistUnsupportedDirective)
{
  SmtSolver s = create_solver(GetParam());
  FunctionalTransitionSystem fts(s);
  EXPECT_THROW(
      SystemVerilogEncoder enc(sv_path("filelist_top.sv"),
                               fts,
                               { sv_path("filelist_bad_directive.f") }),
      PonoException);
}

// A list file naming a source file that does not exist must produce a
// clear parse failure.
TEST_P(SVUnitTests, FilelistMissingFile)
{
  SmtSolver s = create_solver(GetParam());
  FunctionalTransitionSystem fts(s);
  EXPECT_THROW(
      SystemVerilogEncoder enc(
          sv_path("filelist_top.sv"), fts, { sv_path("filelist_missing.f") }),
      PonoException);
}

// Width-parameterized counter.  The property references a
// `localparam` MAX whose value depends on the module parameter,
// exercising NamedValue -> ParameterSymbol resolution.  With
// WIDTH=4 the counter rolls over to MAX=15 at cycle 16 (one reset
// cycle + 15 increments from 0).
TEST_P(SVUnitTests, Parameter) { check_bmc("parameter.sv", 16); }

// Procedural for loop: 4-bit popcount via NamedValue LHS
// accumulation inside an always_comb block.  Falsifies when din
// is all-ones, which lets the registered popcount reach 4 at
// cycle 2 (one reset cycle + one cycle for the always_ff).
TEST_P(SVUnitTests, ForLoop) { check_bmc("for_loop.sv", 2); }

// Generate-for block: four independent counters, each declared
// inside a per-iteration generate block.  The assertion reads
// ctr[2].count via hierarchical reference; the counter reaches 5
// at cycle 6 (one reset cycle + five increments).
TEST_P(SVUnitTests, GenerateBlock) { check_bmc("generate_block.sv", 6); }

// Compound assignment (`|=`) inside a for loop builds an
// OR-reduction.  Exercises ExpressionKind::LValueReference, which
// slang emits as the implicit self-reference on the RHS.  Falsifies
// at cycle 2 once any bit of the free `din` is high.
TEST_P(SVUnitTests, CompoundAssign) { check_bmc("compound_assign.sv", 2); }

// Element-select LHS.  Each bit of a 4-bit register is set
// individually in an unrolled for loop, composing partial writes
// via replace_bits.  Falsifies at cycle 2 when flag reaches
// 4'b1111.
TEST_P(SVUnitTests, ElementSelectLhs) { check_bmc("element_select_lhs.sv", 2); }

// SVA overlapping implication (`|->`): `!rst |-> data != 15`.
// Falsifies at cycle 1 once reset releases and BMC picks data=15.
TEST_P(SVUnitTests, BinaryImplication)
{
  check_bmc("binary_implication.sv", 1);
}

// SVA non-overlapping implication (`|=>`).  The hidden latch is
// 0 at cycle 0 (init), then picks up arm at each cycle.  BMC
// chooses arm=1 during the reset cycle so the latch is high at
// cycle 1 alongside reset_done; with `data` a free input, the
// consequent (data == 10) is falsified at cycle 1.
TEST_P(SVUnitTests, BinaryNonOverlap) { check_bmc("binary_nonoverlap.sv", 1); }

// SVA `and` operator over two simple sub-properties.  Free `data`
// means BMC picks values violating one half immediately at
// cycle 0 (no reset wiring in this design).
TEST_P(SVUnitTests, BinaryAnd) { check_bmc("binary_and.sv", 0); }

// `$past(expr)` returns the previous-cycle value via a hidden
// 1-cycle latch.  reg_data tracks $past(data); the assertion
// `reg_data != $past(data)` is falsified at cycle 1 when BMC
// picks data(0)=0 (so reg_data after reset matches the latched
// past).
TEST_P(SVUnitTests, PastCall) { check_bmc("past_call.sv", 1); }

// `|-> ##N`: the encoder builds N hidden latches to delay the
// antecedent.  With `arm |-> ##2 data == 4'b1010`, BMC picks
// arm=1 during reset (cycle 0); the antecedent latch propagates
// over two cycles and the consequent fails at cycle 2 when BMC
// picks data != 10.
TEST_P(SVUnitTests, SequenceDelay) { check_bmc("sequence_delay.sv", 2); }

// Top-level `always P` unwraps to the same safety property as
// plain `P`.  Counter reaches 5 at cycle 5 (one reset cycle + 5
// increments), falsifying `always (count != 5)`.
TEST_P(SVUnitTests, AlwaysAssertion) { check_bmc("always_assertion.sv", 5); }

// Top-level `s_eventually P` goes through the LTL tableau (its
// negation G!P has no eventuality, so the justice set is just the
// activation latch).  The test harness wraps the TS with
// LivenessToSafetyTranslator and runs BMC.  When `enable` is held low
// after reset, `count` stays at 0 forever and the lasso witness is
// found within a few translator-introduced cycles.
TEST_P(SVUnitTests, EventuallyAssertion)
{
  check_liveness_bmc("eventually_assertion.sv", 5);
}

// Full LTL: request/response liveness `always (req |-> s_eventually
// ack)` nests `always` over `s_eventually`, so it is not pure safety
// and goes through the general tableau.  With `ack` free, BMC finds a
// fair lasso where a request is never acknowledged.
TEST_P(SVUnitTests, ReqAckLiveness) { check_liveness_bmc("req_ack.sv", 10); }

// Same property with `ack` tied high, so it genuinely holds.  The
// tableau must not invent a counterexample; BMC exhausts the bound
// and returns UNKNOWN.
TEST_P(SVUnitTests, ReqAckHolds)
{
  check_liveness_bmc("req_ack_holds.sv", 8, ProverResult::UNKNOWN);
}

// Strong until `busy s_until done` (busy tied high) is a liveness
// obligation that `done` must eventually arrive.  With `done` held
// low forever, the tableau's release-encoded negation G(!done) yields
// a lasso counterexample.
TEST_P(SVUnitTests, UntilLiveness) { check_liveness_bmc("until_live.sv", 10); }

// ---------------------------------------------------------------------------
// Statement kinds
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, IfElse) { check_bmc("if_else.sv", 2); }
TEST_P(SVUnitTests, CaseStatement) { check_bmc("case_stmt.sv", 2); }

// ---------------------------------------------------------------------------
// Expression kinds
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, Ternary) { check_bmc("ternary.sv", 2); }
TEST_P(SVUnitTests, BitSelect) { check_bmc("bit_select.sv", 2); }
TEST_P(SVUnitTests, RangeSelect) { check_bmc("range_select.sv", 2); }
TEST_P(SVUnitTests, Concat) { check_bmc("concat.sv", 2); }
TEST_P(SVUnitTests, Replication) { check_bmc("replication.sv", 2); }

// ---------------------------------------------------------------------------
// Operators
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, BinaryArith) { check_bmc("binary_arith.sv", 2); }
TEST_P(SVUnitTests, BinaryBitwise) { check_bmc("binary_bitwise.sv", 2); }
TEST_P(SVUnitTests, BinaryCompare) { check_bmc("binary_compare.sv", 2); }
TEST_P(SVUnitTests, BinaryLogical) { check_bmc("binary_logical.sv", 2); }
TEST_P(SVUnitTests, Shift) { check_bmc("shift.sv", 2); }
TEST_P(SVUnitTests, UnaryNot) { check_bmc("unary_not.sv", 2); }
TEST_P(SVUnitTests, UnaryMinus) { check_bmc("unary_minus.sv", 2); }
TEST_P(SVUnitTests, Reduction) { check_bmc("reduction.sv", 2); }

// ---------------------------------------------------------------------------
// Combinational logic
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, ContinuousAssign) { check_bmc("continuous_assign.sv", 0); }

TEST_P(SVUnitTests, AlwaysComb) { check_bmc("always_comb.sv", 0); }

TEST_P(SVUnitTests, LegacyAlwaysStar) { check_bmc("legacy_always_comb.sv", 0); }

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVUnitTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests

#endif
