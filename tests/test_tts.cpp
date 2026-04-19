#include <iostream>
#include "core/prop.h"
#include "core/tts.h"
#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "utils/exceptions.h"
#include "engines/bmc.h"
#include "utils/logger.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class TTSUnitTests : public ::testing::Test,
                    public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  void SetUp() override
  {
    s = create_solver(GetParam());
    s2 = create_solver(GetParam());
    tts = new TimedTransitionSystem(s);
    tts2 = new TimedTransitionSystem(s2, TimedAutomatonDelayStrictness::Strict);
    realsort = s->make_sort(smt::REAL);
    intsort = s->make_sort(smt::INT);
    boolsort = s->make_sort(smt::BOOL);
    realsort2 = s2->make_sort(smt::REAL);
    boolsort2 = s2->make_sort(smt::BOOL);
    // logger.set_verbosity(5);
    build_tts();    
  }
  void build_tts() {
    // state variables without state updates
    // will make the system non-functional
    l = tts->make_statevar("l", boolsort);
    b = tts->make_statevar("b", boolsort);
    x = tts->make_statevar("x", realsort);
    y = tts->make_statevar("y", realsort);
    tts->add_nonclock_var(l);
    tts->add_clock_var(x);
    tts->add_clock_var(y);

    // Init state is l /\ b /\ x=0 /\ y=0
    tts->constrain_init(tts->make_term(Equal, x, tts->make_term(0, realsort)));
    tts->constrain_init(tts->make_term(Equal, y, tts->make_term(0, realsort)));
    tts->constrain_init(l);
    tts->constrain_init(b);

    // edge1: (l/\b, x<=2, y:=0, ~l/\b)
    TermVec attributes1 = 
        {b, l, tts->make_term(Le, x, tts->make_term(2, realsort)), // guard x <= 2
        tts->make_term(Equal, tts->next(y), tts->make_term(0, realsort)), // y is reset
        tts->make_term(Equal, tts->next(x), x), // x is unchanged
        tts->make_term(Not, tts->next(l)), // ~l is the next location
        tts->next(b) 
      };
    Term edge1 = s->make_term(And, attributes1);

    // edge2: (~l/\b, y>=2, x:=0, l/\b)
    TermVec attributes2 = 
        {b, tts->make_term(Not, l), tts->make_term(Ge, y, tts->make_term(2, realsort)), // guard y >= 2
        tts->make_term(Equal, tts->next(x), tts->make_term(0, realsort)), // x is reset
        tts->make_term(Equal, tts->next(y), y), // y is unchanged
        tts->next(l), // l is the next location
        tts->next(b) 
      };
    Term edge2 = s->make_term(And, attributes2);
    // edge3: (l/\b, x>=1, y:=0, l/\~b)
    TermVec attributes3 = 
        {b, l, tts->make_term(Ge, x, tts->make_term(1, realsort)),
        tts->make_term(Equal, tts->next(x), x),
        tts->make_term(Equal, tts->next(y), tts->make_term(0, realsort)),
        tts->make_term(Not, tts->next(b)),
        tts->next(l),
      };
    Term edge3 = s->make_term(And, attributes3);

    // edge4: (l/\~b, true, ~l/\~b)
    TermVec attributes4 = 
        {l, tts->make_term(Not, b),
        tts->make_term(Equal, tts->next(x), x),
        tts->make_term(Equal, tts->next(y), y),
        tts->make_term(Not, tts->next(l)),
        tts->make_term(Not, tts->next(b))
      };
    Term edge4 = s->make_term(And, attributes4);

    tts->add_urgent(tts->make_term(And, l, tts->make_term(Not, b)));
    tts->add_invar(
      tts->make_term(Implies, 
        tts->make_term(And, tts->make_term(Not, l), tts->make_term(Not, b)),
        tts->make_term(Le, x, tts->make_term(1, realsort))
      )
    );

    tts->constrain_trans(s->make_term(Or, {edge1, edge2, edge3, edge4}));
    tts->add_invar(s->make_term(Le, s->make_term(Minus, x, y), s->make_term(10, realsort)));
    tts->add_invar(s->make_term(Le, x, s->make_term(Plus, x, s->make_term(10, realsort))));
    tts->add_invar(s->make_term(Le, s->make_term(Minus, x, s->make_term(10, realsort)), x));
    tts->add_invar(s->make_term(Le, x, s->make_term(Minus, s->make_term(10, realsort), s->make_term(0, realsort))));
    tts->add_invar(s->make_term(Implies, l, s->make_term(Le, s->make_term(Minus, x, y), s->make_term(10, realsort))));
    tts->add_invar(s->make_term(Implies, l, s->make_term(Ge, s->make_term(Minus, x, y), s->make_term(-10, realsort))));
    tts->encode_timed_automaton_delays();

    // tts2
    q = tts2->make_statevar("q", boolsort2);
    z = tts2->make_statevar("z", realsort2);
    tts2->add_nonclock_var(q);
    tts2->add_clock_var(z);

    tts2->constrain_init(q);
    tts2->constrain_init(s2->make_term(Equal, z, s2->make_term(0, realsort)));
    Term edge = s2->make_term(And, 
        q, 
        s2->make_term(Not, tts2->next(q)),
        s2->make_term(Equal, z, tts2->next(z))
      );
    Term back_edge = s2->make_term(And, 
        s2->make_term(Not, q), 
        tts2->next(q),
        s2->make_term(Equal, z, tts2->next(z))
      );
    tts2->set_trans(edge);
    tts2->add_invar(s2->make_term(Lt, z, s2->make_term(2, realsort)));
    tts2->add_invar(s2->make_term(Implies, q, s2->make_term(Ge, z, s2->make_term(0, realsort))));
    tts2->encode_timed_automaton_delays();
  }
  SmtSolver s, s2;
  Sort realsort;
  Sort intsort;
  Sort boolsort;
  Sort realsort2;
  Sort boolsort2;
  Term l;
  Term x;
  Term y;
  Term b;
  Term q;
  Term z;
  TimedTransitionSystem * tts;
  TimedTransitionSystem * tts2;
};

TEST_P(TTSUnitTests, TTS_NonEmptyTrans)
{
  EXPECT_FALSE(tts->is_functional());

  Term test_tr1 = s->make_term(And, l, tts->next(s->make_term(Not, l)));
  Term test_tr2 = s->make_term(And, tts->next(l), s->make_term(Not, l));
  
  // Is there a transition from l to ~l, and vice versa?
  s->push();
  s->assert_formula(tts->trans());
  s->assert_formula(test_tr1);
  EXPECT_TRUE(s->check_sat().is_sat());
  s->pop();

  s->push();
  s->assert_formula(tts->trans());
  s->assert_formula(test_tr2);
  EXPECT_TRUE(s->check_sat().is_sat());
  s->pop();

  }
  TEST_P(TTSUnitTests, TTS_BadInvar)
  {
    EXPECT_ANY_THROW(tts->add_invar(x));
    EXPECT_ANY_THROW(tts->add_invar(
        s->make_term(Le, s->make_term(Plus, x, x), s->make_term(0, realsort)))
      );
    EXPECT_ANY_THROW(tts->add_invar(
        s->make_term(And,
          s->make_term(Le, s->make_term(Plus, x, x), s->make_term(0, realsort)),
          l)
        )
      );
  }
  TEST_P(TTSUnitTests, TTS_BMCUnreach)
  {
    Term unreachable_point = 
      tts->make_term(And,
        {tts->make_term(Not, l),
        tts->make_term(Equal, x, tts->make_term(0, realsort)),
        tts->make_term(Equal, y, tts->make_term(3, realsort)),
        b}
      );
    SafetyProperty p1(s, tts->make_term(Not, unreachable_point));
    Bmc b1(p1, *tts, s);
    ProverResult r = b1.check_until(20);
    ASSERT_EQ(r, ProverResult::UNKNOWN);
  }  
  TEST_P(TTSUnitTests, TTS_BMCReach)
  {
    Term reachable_point = 
      tts->make_term(And,
        {b,
        tts->make_term(Not, l),
        tts->make_term(Not,tts->make_term(Equal, x, y)),
        tts->make_term(Ge, x, tts->make_term(1, realsort))
        }
      );
    SafetyProperty p2(s, tts->make_term(Not, reachable_point));
    Bmc b2(p2, *tts, s);
    ProverResult r = b2.check_until(20);
    ASSERT_EQ(r, ProverResult::FALSE);
    std::vector<UnorderedTermMap> cex;
    ASSERT_TRUE(b2.witness(cex));
    ASSERT_TRUE(cex.size() == 3);
  }  

   TEST_P(TTSUnitTests, TTS_BMCReachViaUrgent)
  {
    Term reachable_point = 
      tts->make_term(And,
        {tts->make_term(Not, b),
        tts->make_term(Not, l),
        tts->make_term(Equal,x, tts->make_term(1, realsort))
        }
      );
    SafetyProperty p2(s, tts->make_term(Not, reachable_point));
    Bmc b2(p2, *tts, s);
    ProverResult r = b2.check_until(20);
    ASSERT_EQ(r, ProverResult::FALSE);
    std::vector<UnorderedTermMap> cex;
    ASSERT_TRUE(b2.witness(cex));
    ASSERT_TRUE(cex.size() == 4);
  }

  TEST_P(TTSUnitTests, TTS_BMCUnReachViaUrgent)
  {
    Term unreachable_point = 
      tts->make_term(And,
        {tts->make_term(Not, b),
        l,
        tts->make_term(Gt, y, tts->make_term(0, realsort))
        }
      );
    SafetyProperty p2(s, tts->make_term(Not, unreachable_point));
    Bmc b2(p2, *tts, s);
    ProverResult r = b2.check_until(20);
    ASSERT_EQ(r, ProverResult::UNKNOWN);
  }  

  TEST_P(TTSUnitTests, TTS_BMCUnReachViaInvar)
  {
    Term unreachable_point = 
      tts->make_term(And,
        {tts->make_term(Not, b),
        tts->make_term(Not, l),
        tts->make_term(Gt, x, tts->make_term(1, realsort))
        }
      );
    SafetyProperty p2(s, tts->make_term(Not, unreachable_point));
    Bmc b2(p2, *tts, s);
    ProverResult r = b2.check_until(4);
    ASSERT_EQ(r, ProverResult::UNKNOWN);
  }  


  TEST_P(TTSUnitTests, TTS2_BMCUnReachViaInvar)
  {
    Term unreachable_point = 
      tts2->make_term(And,
        {tts2->make_term(Not, q),
        tts2->make_term(Gt, z, tts2->make_term(2, realsort))
        }
      );
    SafetyProperty p2(s2, tts2->make_term(Not, unreachable_point));
    Bmc b2(p2, *tts2, s2);
    ProverResult r = b2.check_until(3);
    // std::vector<UnorderedTermMap> cex;
    // ASSERT_TRUE(b2.witness(cex));
    // for (auto tv : cex){
    //   for (auto t : tv) {
    //     if (tts->is_curr_var(t.first) || tts->is_input_var(t.first))
    //       std::cout << t.first << " = " << t.second << "\n";
    //   }
    //   std::cout << "\n";
    // }
    ASSERT_EQ(r, ProverResult::UNKNOWN);
  }  

  const std::vector<SolverEnum> real_solver_enums({
      CVC5,
  #ifdef WITH_MSAT
      MSAT,
  #endif
  #ifdef WITH_YICES2
      YICES2,
  #endif
  #ifdef WITH_Z3
      Z3,
  #endif
  });

  INSTANTIATE_TEST_SUITE_P(ParameterizedSolverTTSUnitTests,
                         TTSUnitTests,
                         testing::ValuesIn(real_solver_enums));

}