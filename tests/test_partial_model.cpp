#include "gtest/gtest.h"
#include "utils/partial_model.h"
#include "utils/sygus_ic3formula_helper.h"
#include "smt/available_solvers.h"


using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class DynamicCoiUnitTests : 
                     public ::testing::Test,
                     public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  void SetUp() override
  {
    s = create_solver(GetParam());
    s->set_opt("produce-models", "true");
    s->set_opt("incremental", "true");
    boolsort = s->make_sort(BOOL);
    bvsort8 = s->make_sort(BV, 8);
  }
  SmtSolver s;
  Sort boolsort, bvsort8;
};


#define NOT(x)    (s->make_term(Not, (x)))
#define ADD(x, y) (s->make_term(BVAdd, (x), (y)))
#define SUB(x, y) (s->make_term(BVSub, (x), (y)))
#define BVAND(x, y) (s->make_term(BVAnd, (x), (y)))
#define BVOR(x, y) (s->make_term(BVOr, (x), (y)))
#define EQ(x, y) (s->make_term(Equal, (x), (y)))
//#define BoolEQ(x, y) (s->make_term(Iff, (x), (y)))
#define ITE(c, x, y) (s->make_term(Ite, (c), (x), (y)))

#define CheckPartialModel(p,u) \
    {      \
      s->push();                         \
      auto ast = EQ(p, u);        \
      s->assert_formula(ast);            \
      if ( s->check_sat().is_sat() ) {   \
        auto m_cube = pt.GetPartialModelInCube(ast);   \
        std::cout << "expr: " << ast << std::endl; \
        std::cout << m_cube.first.term << std::endl;     \
        std::cout << m_cube.second.to_string() << std::endl; \
        s->pop(); s->push();\
        s->assert_formula(NOT(ast));\
        s->assert_formula(m_cube.first.term);\
        EXPECT_TRUE(s->check_sat().is_unsat()); \
      }                                  \
      s->pop();                          \
    }

TEST_P(DynamicCoiUnitTests, SimpleCoiTest)
{
  
    Term a = s->make_symbol("a", bvsort8);
    Term b = s->make_symbol("b", bvsort8);
    Term c = s->make_symbol("c", bvsort8);
    Term d = s->make_symbol("d", bvsort8);
    Term e = s->make_symbol("e", bvsort8);
    Term f = s->make_symbol("f", bvsort8);
    Term g = s->make_symbol("g", bvsort8);

    Term u = s->make_symbol("u", bvsort8);
    Term v = s->make_symbol("v", bvsort8);
    Term w = s->make_symbol("w", bvsort8);
    Term x = s->make_symbol("x", bvsort8);
    Term y = s->make_symbol("y", bvsort8);
    Term z = s->make_symbol("z", bvsort8);
    
    auto x_plus_1 = s->make_term(BVAdd, x, s->make_term(1, bvsort8));
    auto x_sub_1 = s->make_term(BVSub, x, s->make_term(1, bvsort8));
    auto x_and_1 = s->make_term(BVAnd, x, s->make_term(1, bvsort8));
    auto x_plus_y = ADD(x, y);
    auto x_sub_y = SUB(x, y);
    auto x_and_y = BVAND(x, y);
    auto t0 = BVAND(x, x_sub_1);
    auto t1 = BVAND(y, x_sub_1);
    auto t2 = ADD(x, x_sub_1);
    auto t3 = ADD(x, s->make_term(BVSub, x, s->make_term(1, bvsort8)));
    auto t4 = ADD(x, x_and_1);
    auto t5 = ADD(y, x_and_1);
    auto t6 = ADD(y, s->make_term(BVAnd, x, s->make_term(1, bvsort8)));

    auto e1 = ITE(EQ(SUB(x, y), z), ADD(a, b), c);
    auto e2 = ITE(EQ(SUB(x, y), x), ADD(a, b), SUB(a, c));
    auto e3 = ITE(EQ(BVAND(x, y), x), BVAND(a, b), BVOR(a, c));
    auto e4a = ITE(EQ(SUB(u, SUB(v,w)), d), ADD(e, f), g);
    auto e4b = ITE( EQ( EQ(BVAND(a, SUB(e1,e2)), f), EQ(BVAND(x, y), x) ), ADD(e4a, f), e1);


    PartialModelGen pt(s);
    
    CheckPartialModel(x_plus_1, u);
    CheckPartialModel(x_sub_1, u);
    CheckPartialModel(x_and_1, u);
    CheckPartialModel(x_plus_y, u);
    CheckPartialModel(x_sub_y, u);
    CheckPartialModel(x_and_y, u);
    CheckPartialModel(t0, u);
    CheckPartialModel(t1, u);
    CheckPartialModel(t2, u);
    CheckPartialModel(t3, u);
    CheckPartialModel(t4, u);
    CheckPartialModel(t5, u);
    CheckPartialModel(t6, u);
    CheckPartialModel(e1, u);
    CheckPartialModel(e2, u);
    CheckPartialModel(e3, u);
    CheckPartialModel(e4a, u);
    CheckPartialModel(e4b, u);
}


INSTANTIATE_TEST_SUITE_P(ParameterizedDynamicCoiUnitTests,
                         DynamicCoiUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
