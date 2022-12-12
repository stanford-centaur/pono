#include <string>
#include <tuple>
#include <vector>

#include "core/fts.h"
#include "frontends/btor2_encoder.h"
#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "test_encoder_inputs.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

// BTOR2 systems that are too large to check semantic equality
// of the terms
const unordered_set<string> large_files({ "ridecore.btor" });

// Test that copying TransitionSystems to different solvers
// doesn't change the semantics of the TransitionSystem
class CopyUnitTests
    : public ::testing::Test,
      public ::testing::WithParamInterface<tuple<SolverEnum, string>>
{
};

TEST_P(CopyUnitTests, CopyFromDefault)
{
  cout << "Copy to " << get<0>(GetParam()) << endl;
  FunctionalTransitionSystem fts;
  fts.solver()->set_opt("incremental", "true");
  // PONO_SRC_DIR is a macro set using CMake PROJECT_SRC_DIR
  string filename = STRFY(PONO_SRC_DIR);
  filename += "/tests/encoders/inputs/btor2/";
  filename += get<1>(GetParam());
  cout << "Reading file: " << filename << ", for copy test" << endl;
  BTOR2Encoder be(filename, fts);

  SmtSolver s = create_solver(get<0>(GetParam()));
  // a term translator that can translate from one other solver to this solver
  TermTranslator tt(s);

  // try copying the transition system
  FunctionalTransitionSystem fts_copy(fts, tt);

  EXPECT_EQ(fts.statevars().size(), fts_copy.statevars().size());
  EXPECT_EQ(fts.inputvars().size(), fts_copy.inputvars().size());
  EXPECT_EQ(fts.named_terms().size(), fts_copy.named_terms().size());
  EXPECT_EQ(fts.state_updates().size(), fts_copy.state_updates().size());

  // now copy it back
  // first need a TermTranslator with the cache populated for the symbols
  TermTranslator back_tt(fts.solver());
  UnorderedTermMap & cache = back_tt.get_cache();
  for (auto v : fts.statevars()) {
    cache[tt.transfer_term(v)] = v;
    cache[tt.transfer_term(fts.next(v))] = fts.next(v);
  }
  for (auto v : fts.inputvars()) {
    cache[tt.transfer_term(v)] = v;
    cache[tt.transfer_term(fts.next(v))] = fts.next(v);
  }

  FunctionalTransitionSystem fts_2(fts_copy, back_tt);
  ASSERT_EQ(fts.is_functional(), fts_2.is_functional());
  const SmtSolver & fts_solver = fts.solver();
  ASSERT_EQ(fts_solver, fts_2.solver());

  // VARS -- these should be exactly the same
  for (auto v : fts.inputvars()) {
    EXPECT_TRUE(fts_2.inputvars().find(v) != fts_2.inputvars().end());
  }
  for (auto v : fts.statevars()) {
    EXPECT_TRUE(fts_2.statevars().find(v) != fts_2.statevars().end());
    EXPECT_TRUE(fts.next(v) == fts_2.next(v));
  }

  // now check that each of the data structures were copied to the
  // TransitionSystem correctly they could be syntactically different, but they
  // should be semantically equivalent
  Result r;

  // INIT
  fts_solver->push();
  fts_solver->assert_formula(
      fts_solver->make_term(Distinct, fts.init(), fts_2.init()));
  r = fts_solver->check_sat();
  EXPECT_TRUE(r.is_unsat());
  fts_solver->pop();

  // this check is too expensive on large systems
  if (large_files.find(get<1>(GetParam())) == large_files.end()) {
    // TRANS
    fts_solver->push();
    fts_solver->assert_formula(
        fts_solver->make_term(Distinct, fts.trans(), fts_2.trans()));
    r = fts_solver->check_sat();
    EXPECT_TRUE(r.is_unsat());
    fts_solver->pop();
  }
}

TEST_P(CopyUnitTests, CopyToDefault)
{
  cout << "Copy from " << get<0>(GetParam()) << endl;
  SmtSolver s = create_solver(get<0>(GetParam()));
  s->set_opt("incremental", "true");
  FunctionalTransitionSystem fts(s);
  // PONO_SRC_DIR is a macro set using CMake PROJECT_SRC_DIR
  string filename = STRFY(PONO_SRC_DIR);
  filename += "/tests/encoders/inputs/btor2/";
  filename += get<1>(GetParam());
  cout << "Reading file: " << filename << ", for copy test" << endl;
  BTOR2Encoder be(filename, fts);

  // try copying the transition system
  FunctionalTransitionSystem fts_copy;

  // a term translator that can translate from one other solver to this solver
  TermTranslator tt(fts_copy.solver());

  fts_copy = FunctionalTransitionSystem(fts, tt);

  EXPECT_EQ(fts.statevars().size(), fts_copy.statevars().size());
  EXPECT_EQ(fts.inputvars().size(), fts_copy.inputvars().size());
  EXPECT_EQ(fts.named_terms().size(), fts_copy.named_terms().size());
  EXPECT_EQ(fts.state_updates().size(), fts_copy.state_updates().size());

  // now copy it back
  // first need a TermTranslator with the cache populated for the symbols
  ASSERT_EQ(fts.solver(), s);
  TermTranslator back_tt(s);
  UnorderedTermMap & cache = back_tt.get_cache();
  for (auto v : fts.statevars()) {
    cache[tt.transfer_term(v)] = v;
    cache[tt.transfer_term(fts.next(v))] = fts.next(v);
  }
  for (auto v : fts.inputvars()) {
    cache[tt.transfer_term(v)] = v;
    cache[tt.transfer_term(fts.next(v))] = fts.next(v);
  }

  FunctionalTransitionSystem fts_2(fts_copy, back_tt);
  ASSERT_EQ(fts.is_functional(), fts_2.is_functional());
  ASSERT_EQ(s, fts_2.solver());

  // VARS -- these should be exactly the same
  for (auto v : fts.inputvars()) {
    EXPECT_TRUE(fts_2.inputvars().find(v) != fts_2.inputvars().end());
  }
  for (auto v : fts.statevars()) {
    EXPECT_TRUE(fts_2.statevars().find(v) != fts_2.statevars().end());
    EXPECT_TRUE(fts.next(v) == fts_2.next(v));
  }

  // now check that each of the data structures were copied to the
  // TransitionSystem correctly they could be syntactically different, but they
  // should be semantically equivalent
  Result r;

  // INIT
  s->push();
  s->assert_formula(s->make_term(Distinct, fts.init(), fts_2.init()));
  r = s->check_sat();
  EXPECT_TRUE(r.is_unsat());
  s->pop();

  // HACK -- only check large systems with boolector
  // otherwise the tests take too long for large files
  if (large_files.find(get<1>(GetParam())) == large_files.end()
      || s->get_solver_enum() == smt::BTOR) {
    // TRANS
    s->push();
    s->assert_formula(s->make_term(Distinct, fts.trans(), fts_2.trans()));
    r = s->check_sat();
    EXPECT_TRUE(r.is_unsat());
    s->pop();
  }
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedSolverCopyUnitTests,
    CopyUnitTests,
    testing::Combine(
        testing::ValuesIn(available_solver_enums()),
        // from test_encoder_inputs.h
        testing::ValuesIn(btor2_inputs)));

}  // namespace pono_tests
