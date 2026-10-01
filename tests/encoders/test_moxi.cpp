#include <cstdio>
#include <fstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "core/rts.h"
#include "engines/kinduction.h"
#include "frontends/moxi_encoder.h"
#include "gtest/gtest.h"
#include "modifiers/prop_monitor.h"
#include "smt/available_solvers.h"
#include "test_encoder_inputs.h"
#include "utils/exceptions.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

string moxi_path(const string & file)
{
  // PONO_SRC_DIR is a macro set using CMake PROJECT_SRC_DIR
  return string(STRFY(PONO_SRC_DIR)) + "/tests/encoders/inputs/moxi/" + file;
}

/** @return whether the two formulas are equivalent */
bool equivalent(const SmtSolver & s, const Term & a, const Term & b)
{
  s->push();
  s->assert_formula(s->make_term(Distinct, a, b));
  const bool result = s->check_sat().is_unsat();
  s->pop();
  return result;
}

/** Each query with each solver that has the theories it needs. */
vector<tuple<SolverEnum, MoxiQuery>> moxi_query_params()
{
  vector<tuple<SolverEnum, MoxiQuery>> params;
  for (const MoxiQuery & query : moxi_queries) {
    for (const SolverEnum se : filter_solver_enums(query.theories)) {
      if (!query.excluded.count(se)) {
        params.emplace_back(se, query);
      }
    }
  }
  return params;
}

class MoxiQueryUnitTests
    : public ::testing::Test,
      public ::testing::WithParamInterface<tuple<SolverEnum, MoxiQuery>>
{
};

TEST_P(MoxiQueryUnitTests, Encode)
{
  const SolverEnum se = get<0>(GetParam());
  const MoxiQuery & query = get<1>(GetParam());
  SmtSolver s = create_solver(se);
  RelationalTransitionSystem rts(s);
  MoxiEncoder encoder(moxi_path(query.file), rts, query.index);

  Term prop = encoder.prop();
  // A reachability condition with next-state variables is monitored, as
  // pono does it too.
  if (!rts.only_curr(prop)) {
    prop = add_prop_monitor(rts, prop);
  }
  SafetyProperty p(s, prop);
  KInduction kind(p, rts, s);
  EXPECT_EQ(kind.check_until(20), query.result);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverMoxiQueryUnitTests,
                         MoxiQueryUnitTests,
                         testing::ValuesIn(moxi_query_params()));

class MoxiEncoderUnitTests : public ::testing::Test
{
 protected:
  // cvc5 comes with every build and has all the theories the inputs use.
  void SetUp() override { s = create_solver(CVC5); }

  SmtSolver s;
};

TEST_F(MoxiEncoderUnitTests, SelectsACheckSystemCommand)
{
  RelationalTransitionSystem rts(s);
  MoxiEncoder encoder(moxi_path("checks.moxi"), rts, 1);
  EXPECT_EQ(encoder.query_name(), "impossible");
}

// Only the command with the :queries attribute is rejected; moxi_queries
// checks the others of the file.
TEST_F(MoxiEncoderUnitTests, RejectsTheQueriesAttributeOfTheSelectedCommand)
{
  RelationalTransitionSystem rts(s);
  try {
    MoxiEncoder encoder(moxi_path("checks.moxi"), rts, 3);
    FAIL() << "the :queries attribute was encoded";
  }
  catch (const PonoException & e) {
    EXPECT_NE(string(e.what()).find(
                  "checks.moxi:31:3: the :queries attribute is not supported"),
              string::npos)
        << e.what();
  }
}

TEST_F(MoxiEncoderUnitTests, NamesTheVariablesOfSubsystemsAfterTheirInstance)
{
  RelationalTransitionSystem rts(s);
  MoxiEncoder encoder(moxi_path("subsystems.moxi"), rts);
  // The check-system renames the variables of the system it checks.
  for (const string name : { "x",
                             "y",
                             "m",
                             "steps",
                             "d1.mid",
                             "d2.mid",
                             "d1.first.held",
                             "d2.second.held" }) {
    EXPECT_NO_THROW(rts.lookup(name)) << name;
  }
  for (const string name : { "in", "out", "mid", "t", "held" }) {
    EXPECT_THROW(rts.lookup(name), PonoException) << name;
  }
}

// A variable that only the transitions refer to, unprimed, need not be a
// state, which depends on the query as its assumptions refer to variables too.
TEST_F(MoxiEncoderUnitTests, MakesInputsOfVariablesOnlyTransitionsReferTo)
{
  RelationalTransitionSystem free_rts(s);
  MoxiEncoder free(moxi_path("assumptions.moxi"), free_rts, 0);
  EXPECT_TRUE(free_rts.is_input_var(free_rts.lookup("i")));
  EXPECT_TRUE(free_rts.is_curr_var(free_rts.lookup("s")));

  // a solver of its own, as it cannot have two variables of the same name
  RelationalTransitionSystem assumed_rts(create_solver(CVC5));
  MoxiEncoder assumed(moxi_path("assumptions.moxi"), assumed_rts, 1);
  EXPECT_TRUE(assumed_rts.is_curr_var(assumed_rts.lookup("i")));
}

// Primed, an input refers to its next value, so it has to be a state.
TEST_F(MoxiEncoderUnitTests, MakesStatesOfPrimedInputs)
{
  RelationalTransitionSystem rts(s);
  MoxiEncoder encoder(moxi_path("assumptions.moxi"), rts, 3);
  EXPECT_TRUE(rts.is_curr_var(rts.lookup("b")));
}

TEST_F(MoxiEncoderUnitTests, KeepsDeclaredConstants)
{
  RelationalTransitionSystem rts(s);
  MoxiEncoder encoder(moxi_path("constants.moxi"), rts);
  const Term k = rts.lookup("k");
  ASSERT_TRUE(rts.is_curr_var(k));
  EXPECT_EQ(rts.state_updates().at(k), k);
}

// A :current formula that the query lists replaces the initial conditions of
// the system and of its subsystems, rather than joining them.
TEST_F(MoxiEncoderUnitTests, ReplacesTheInitialConditionsWithTheCurrentFormula)
{
  const string path = testing::TempDir() + "pono_current.moxi";
  {
    ofstream file(path);
    file << "(set-logic QF_LIA)\n"
         << "(define-system Inner :output ((y Int)) :init (= y 1)"
         << " :trans (= y' y))\n"
         << "(define-system Outer :output ((x Int)) :local ((y Int))"
         << " :init (= x 0) :trans (= x' (+ x y)) :subsys (sub (Inner y)))\n";
    // the same formulas twice, of which only the first query lists :current
    for (const string formulas : { "(high seven)", "(seven)" }) {
      file << "(check-system Outer :output ((x Int)) :local ((y Int))"
           << " :current (high (> x 5)) :reachable (seven (= x 7))"
           << " :query (q " << formulas << "))\n";
    }
  }
  RelationalTransitionSystem current_rts(s);
  MoxiEncoder current(path, current_rts, 0);
  // a solver of its own, as it cannot have two variables of the same name
  SmtSolver s2 = create_solver(CVC5);
  RelationalTransitionSystem init_rts(s2);
  MoxiEncoder init(path, init_rts, 1);
  remove(path.c_str());

  EXPECT_TRUE(equivalent(
      s,
      current_rts.init(),
      s->make_term(
          Gt, current_rts.lookup("x"), s->make_term(5, s->make_sort(INT)))));
  const Sort int_sort = s2->make_sort(INT);
  EXPECT_TRUE(equivalent(
      s2,
      init_rts.init(),
      s2->make_term(
          And,
          s2->make_term(
              Equal, init_rts.lookup("x"), s2->make_term(0, int_sort)),
          s2->make_term(
              Equal, init_rts.lookup("y"), s2->make_term(1, int_sort)))));
}

// The symbols standing for variables while a file is read are apart from
// those of any other file.
TEST_F(MoxiEncoderUnitTests, SharesASolverBetweenFiles)
{
  RelationalTransitionSystem counter_rts(s);
  MoxiEncoder counter(moxi_path("counter.moxi"), counter_rts);
  RelationalTransitionSystem current_rts(s);
  EXPECT_NO_THROW(MoxiEncoder(moxi_path("current.moxi"), current_rts));
}

TEST_F(MoxiEncoderUnitTests, RejectsACheckSystemIndexOutOfRange)
{
  RelationalTransitionSystem rts(s);
  try {
    MoxiEncoder encoder(moxi_path("counter.moxi"), rts, 4);
    FAIL() << "check-system 4 of 4 was encoded";
  }
  catch (const PonoException & e) {
    EXPECT_NE(string(e.what()).find("out of range"), string::npos) << e.what();
  }
}

/** A malformed MoXI file, with part of the message rejecting it. */
struct MoxiError
{
  string text;
  string message;
};

ostream & operator<<(ostream & os, const MoxiError & error)
{
  return os << error.message;
}

// Each file is well-formed apart from one mistake, mostly the one that the
// reference sort checker in the MoXI tool suite tests with the same name.
const vector<MoxiError> moxi_errors({
    // lexical rules
    { "(define-fun f () Int 00123)", "numeral with a leading zero" },
    { "(set-logic QF_LRA) (define-fun f () Real 0.)", "syntax error" },
    { "(define-system S :input (('abc Bool)))", "unexpected character" },
    { "(define-system S :input ((a'bc Bool)))", "syntax error" },
    { "(define-fun |abc () Bool true)", "without a closing |" },
    { "(define-fun |abc\\| () Bool true)", "without a closing |" },
    // commands
    { "(fake-command)", "unknown or unsupported command fake-command" },
    { "(set-logic PHONY)", "unknown logic PHONY" },
    { "(set-logic QF_BV) (set-logic QF_BV)", "the logic is already set" },
    { "(define-fun f () Bool true) (define-fun f () Bool false)",
      "symbol f is already declared" },
    { "(declare-enum-sort S (and))", "and is a predefined symbol" },
    { "(declare-enum-sort Bool (A B C))", "sort Bool is already defined" },
    { "(set-logic QF_LIA) (define-fun f ((x Int) (y Int)) Int (> x y))",
      "the definition of f has sort Bool, but Int is declared" },
    { "(declare-fun x (Bool Bool) Bool) (define-fun y () Bool (x false))",
      "function x takes 2 arguments, not 1" },
    { "(set-logic QF_BV) (define-fun f () (_ BitVec 8) (bvadd #x01 #b1))",
      "cannot apply bvadd" },
    { "(define-fun f () Bool (to_bv3 true))", "unknown function to_bv3" },
    { "(define-fun f () Bool g)", "unknown symbol g" },
    { "(define-fun f () Bool (let ((y true)) y'))",
      "the bound variable y cannot be primed" },
    { "(define-fun c () Bool true) (define-system S :init c')",
      "only variables of a system can be primed" },
    // define-system
    { "(define-system S :output ((o Bool)) :init o')",
      "the next-state variable o' cannot appear in :init" },
    { "(define-system S :output ((o Bool)) :inv o')",
      "the next-state variable o' cannot appear in :inv" },
    { "(define-system S :input ((i Bool)) :output ((i Bool)))",
      "variable i is declared more than once" },
    { "(define-system S :input ((i Bool)) :input ((j Bool)))",
      "attribute :input appears more than once" },
    { "(define-system S :output ((o Bool)) :init o :init (not o))",
      "attribute :init appears more than once" },
    { "(define-system S :init true :input ((i Bool)))", "syntax error" },
    { "(define-system S :init (bvnot #b0))", "instead of Bool" },
    { "(define-system S1 :input ((i Bool)))"
      "(define-system S2 :input ((i Bool)) :subsys (A (S0 i)))",
      "unknown system S0" },
    { "(define-system S1 :input ((i Bool)) :output ((o Bool)))"
      "(define-system S2 :input ((i Bool)) :subsys (A (S1 i)))",
      "passes 1 variables to system S1, which has 2 inputs and outputs" },
    { "(define-system S1 :input ((i Bool)))"
      "(define-system S2 :subsys (A (S1 in)))",
      "passes in, which is not a variable of system S2" },
    { "(set-logic QF_BV) (define-system S1 :input ((i Bool)))"
      "(define-system S2 :local ((b (_ BitVec 1))) :subsys (A (S1 b)))",
      "passes b of sort (_ BitVec 1) for i of sort Bool" },
    { "(define-system S1) (define-system S2 :local ((l Bool)) :subsys (l "
      "(S1)))",
      "subsystem l has the name of another variable or subsystem of S2" },
    // check-system
    { "(define-system S) (check-system S0)", "unknown system S0" },
    { "(define-system S :input ((i Bool))) (check-system S :input ())",
      "the :input attribute declares 0 variables, but system S has 1" },
    { "(set-logic QF_LIA) (define-system S :input ((i Bool)))"
      "(check-system S :input ((i Int)))",
      "variable i has sort Int, but the variable i of system S it renames "
      "has sort Bool" },
    { "(define-system S :output ((o Bool))) (check-system S :current (c o'))",
      "the next-state variable o' cannot appear in :current" },
    { "(define-system S :output ((o Bool)))"
      "(check-system S :assumption (a o) :fairness (a o))",
      "another formula of this check-system is named a" },
    { "(define-system S :output ((o Bool)))"
      "(check-system S :reachable (r o) :query (q (f)))",
      "query q refers to f, which is not a formula of this check-system" },
    { "(define-system S) (check-system S :query (q ()))",
      "query q lists no formulas" },
    { "(define-system S) (check-system S :queries ())", "syntax error" },
    { "(define-system S :output ((o Bool)))"
      "(check-system S :reachable (r o) :query (q (r)) :query (q (r)))",
      "query q is defined more than once" },
    { "(define-system S :output ((o Bool)))"
      "(check-system S :current (c1 o) :current (c2 (not o))"
      " :reachable (r o) :query (q (c1 c2 r)))",
      "query q has more than one :current formula" },
    // what the encoder does not support
    { "(define-system S :output ((o Bool)))"
      "(check-system S :fairness (f o) :reachable (r o) :query (q (f r)))",
      "only queries without fairness conditions are supported" },
    { "(define-system S :output ((o Bool)))"
      "(check-system S :reachable (r o) :queries ((q1 (r)) (q2 (r))))",
      "the :queries attribute is not supported" },
    { "(define-system S :output ((o Bool)))", "has no check-system command" },
    { "(define-system S) (check-system S)",
      "the check-system command has no query" },
    { "(define-system S :output ((o Bool)))"
      "(check-system S :reachable (r o) :query (q1 (r)) :query (q2 (r)))",
      "has 2 queries, but only one query per command is supported" },
});

class MoxiErrorUnitTests : public ::testing::Test,
                           public ::testing::WithParamInterface<MoxiError>
{
};

TEST_P(MoxiErrorUnitTests, Rejects)
{
  const string path = testing::TempDir() + "pono_test.moxi";
  {
    ofstream file(path);
    file << GetParam().text << "\n";
  }
  SmtSolver s = create_solver(CVC5);
  RelationalTransitionSystem rts(s);
  try {
    MoxiEncoder encoder(path, rts);
    ADD_FAILURE() << "the file was encoded:\n" << GetParam().text;
  }
  catch (const PonoException & e) {
    EXPECT_NE(string(e.what()).find(GetParam().message), string::npos)
        << e.what();
  }
  remove(path.c_str());
}

INSTANTIATE_TEST_SUITE_P(MoxiErrorUnitTests,
                         MoxiErrorUnitTests,
                         testing::ValuesIn(moxi_errors));

// An error names the file and the line and column where it occurs.
TEST(MoxiLocationUnitTests, LocatesTheError)
{
  const string path = testing::TempDir() + "pono_test.moxi";
  {
    ofstream file(path);
    file << "(set-logic QF_BV)\n"
         << "; a comment\n"
         << "(define-system S\n"
         << "  :output ((o Bool))\n"
         << "  :init (and o p))\n";
  }
  SmtSolver s = create_solver(CVC5);
  RelationalTransitionSystem rts(s);
  try {
    MoxiEncoder encoder(path, rts);
    FAIL() << "the file was encoded";
  }
  catch (const PonoException & e) {
    EXPECT_EQ(string(e.what()), path + ":5:16: unknown symbol p");
  }
  remove(path.c_str());
}

}  // namespace pono_tests
