#pragma once

#ifdef WITH_SLANG

#include <string>
#include <vector>

#include "core/fts.h"
#include "engines/bmc.h"
#include "frontends/systemverilog/encoder.h"
#include "gtest/gtest.h"
#include "modifiers/control_signals.h"
#include "modifiers/liveness_to_safety_translator.h"
#include "modifiers/mod_ts_prop.h"
#include "smt-switch/utils.h"
#include "smt/available_solvers.h"
#include "test_encoder_inputs.h"
#include "utils/exceptions.h"

// used in the individual tests to recover string paths from macros
#ifndef STRHELPER
#define STRHELPER(A) #A
#define STRFY(A) STRHELPER(A)
#endif

namespace pono_tests {

// Shared fixture for every SystemVerilogEncoder test file.  Split out of the
// original test_systemverilog.cpp so the thematic test files
// (test_systemverilog_{types,operators,statements,hierarchy,sva,unsupported,
// integration}.cpp) can each `#include` it and run their own
// INSTANTIATE_TEST_SUITE_P against the same fixture, instead of duplicating
// sv_path()/find_reset()/check_bmc().
class SVUnitTests : public ::testing::Test,
                    public ::testing::WithParamInterface<smt::SolverEnum>
{
 protected:
  static std::string sv_path(const std::string & name)
  {
    return std::string(STRFY(PONO_SRC_DIR))
           + "/tests/encoders/inputs/systemverilog/" + name;
  }

  // Locate the `rst` input port in the encoded transition system,
  // returning a null Term if the design has no such input.  The
  // encoder names ports as "<top_module>.<port>", so we scan
  // named_terms for any entry ending in ".rst".
  static smt::Term find_reset(const pono::TransitionSystem & ts)
  {
    const std::string suffix = ".rst";
    for (const auto & [name, term] : ts.named_terms()) {
      if (name.size() >= suffix.size()
          && name.compare(name.size() - suffix.size(), suffix.size(), suffix)
                 == 0) {
        return term;
      }
    }
    return smt::Term();
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
  void check_bmc(const std::string & file,
                 size_t bound,
                 pono::ProverResult expected = pono::ProverResult::FALSE,
                 const std::vector<std::string> & filelists = {})
  {
    check_prover<pono::Bmc>(file, bound, expected, filelists);
  }

  // Generalization of check_bmc() over any SafetyProver-shaped engine
  // (Bmc, KInduction, IC3Bits, ...) that shares the EngineT(prop, ts,
  // solver) constructor and the check_until()/witness_length()
  // interface.  Used by the integration tests to exercise more than
  // just Bmc against the SV frontend's output.
  template <typename EngineT>
  void check_prover(const std::string & file,
                    size_t bound,
                    pono::ProverResult expected = pono::ProverResult::FALSE,
                    const std::vector<std::string> & filelists = {})
  {
    using namespace pono;
    using namespace smt;

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
    EngineT engine(prop, ts, s);
    EXPECT_EQ(engine.check_until(bound), expected);
    if (expected == ProverResult::FALSE) {
      EXPECT_EQ(engine.witness_length(), bound);
    }
  }

  // Encode `file` (which must contain a single temporal/LTL
  // assertion), reduce its generalized-Büchi justice set to a safety
  // property via LivenessToSafetyTranslator, then run BMC up to
  // `bound`.  Reset preprocessing is applied if the design has an
  // `rst` input.  The L2S translator adds a "save" input and a few
  // bookkeeping state vars that BMC chooses freely; a FALSE verdict
  // means BMC found a fair lasso violating the LTL property.
  void check_liveness_bmc(
      const std::string & file,
      size_t bound,
      pono::ProverResult expected = pono::ProverResult::FALSE)
  {
    using namespace pono;
    using namespace smt;

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

  // Assert that encoding `file` throws a PonoException (the clean-
  // rejection contract used by type_to_sort() and the
  // BinaryOperator/UnaryOperator default cases for genuinely
  // out-of-scope constructs).
  void expect_encode_throws(const std::string & file,
                            const std::vector<std::string> & filelists = {})
  {
    using namespace pono;
    using namespace smt;
    SmtSolver s = create_solver(GetParam());
    FunctionalTransitionSystem fts(s);
    EXPECT_THROW(SystemVerilogEncoder enc(sv_path(file), fts, filelists),
                 PonoException);
  }

  // Assert that encoding `file` succeeds without throwing at all --
  // used for out-of-scope constructs that turn out (verified
  // empirically) not to hit any of the encoder's rejection paths:
  // the construct's containing symbol/statement is never walked in
  // the first place (e.g. `program`/`checker` instances, `bind`,
  // `defparam`, `fork`/`wait` statements), so encoding quietly
  // proceeds as if the construct weren't there.  Distinct from
  // expect_silently_dropped(), which is about one *assertion*
  // vanishing from propvec()/ltl_justice() -- here the whole
  // construct (often not an assertion at all) has no observable
  // effect on the resulting transition system whatsoever.
  void expect_encode_succeeds_ignoring(
      const std::string & file, const std::vector<std::string> & filelists = {})
  {
    using namespace pono;
    using namespace smt;
    SmtSolver s = create_solver(GetParam());
    FunctionalTransitionSystem fts(s);
    EXPECT_NO_THROW(SystemVerilogEncoder enc(sv_path(file), fts, filelists));
  }
};

}  // namespace pono_tests

#endif  // WITH_SLANG
