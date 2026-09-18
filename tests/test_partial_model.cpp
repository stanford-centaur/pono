#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "utils/partial_model.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class DynamicCoiUnitTests : public ::testing::Test,
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

#define NOT(x) (s->make_term(Not, (x)))
#define ADD(x, y) (s->make_term(BVAdd, (x), (y)))
#define SUB(x, y) (s->make_term(BVSub, (x), (y)))
#define BVAND(x, y) (s->make_term(BVAnd, (x), (y)))
#define BVOR(x, y) (s->make_term(BVOr, (x), (y)))
#define EQ(x, y) (s->make_term(Equal, (x), (y)))
#define BVULT(x, y) (s->make_term(BVUlt, (x), (y)))
#define L_AND(x, y) (s->make_term(And, (x), (y)))
// #define BoolEQ(x, y) (s->make_term(Iff, (x), (y)))
#define ITE(c, x, y) (s->make_term(Ite, (c), (x), (y)))

#define CheckPartialModel(p, u)                    \
  {                                                \
    s->push();                                     \
    auto ast = EQ(p, u);                           \
    s->assert_formula(ast);                        \
    if (s->check_sat().is_sat()) {                 \
      auto m_cube = pt.GetPartialModelInCube(ast); \
      s->pop();                                    \
      s->push();                                   \
      s->assert_formula(NOT(ast));                 \
      s->assert_formula(m_cube.first.term);        \
      EXPECT_TRUE(s->check_sat().is_unsat());      \
    }                                              \
    s->pop();                                      \
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
  auto e4a = ITE(EQ(SUB(u, SUB(v, w)), d), ADD(e, f), g);
  auto e4b = ITE(
      EQ(EQ(BVAND(a, SUB(e1, e2)), f), EQ(BVAND(x, y), x)), ADD(e4a, f), e1);

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

#define CheckPartialModel_bitlevel(ast_in)              \
  {                                                     \
    s->push();                                          \
    auto ast = (ast_in);                                \
    s->assert_formula(ast);                             \
    if (s->check_sat().is_sat()) {                      \
      auto m_cube = pt.GetPartialModelInCube_bitlevel(ast);  \
      s->pop();                                         \
      s->push();                                        \
      s->assert_formula(NOT(ast));                      \
      s->assert_formula(m_cube.first.term);             \
      EXPECT_TRUE(s->check_sat().is_unsat());           \
    }                                                   \
    s->pop();                                           \
  }

TEST_P(DynamicCoiUnitTests, BitCoiTest)
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

  auto c0 = s->make_term(0, bvsort8);
  auto c1 = s->make_term(1, bvsort8);
  auto c2 = s->make_term(2, bvsort8);
  auto c4 = s->make_term(4, bvsort8);
  auto c8 = s->make_term(8, bvsort8);

  auto x_plus_1 = ADD(x, c1);
  auto x_sub_1 = SUB(x, c1);
  auto x_and_1 = BVAND(x, c1);
  auto x_plus_y = ADD(x, y);
  auto x_sub_y = SUB(x, y);
  auto x_and_y = BVAND(x, y);
  auto t0 = BVAND(x, x_sub_1);
  auto t1 = BVAND(y, x_sub_1);
  auto t2 = ADD(x, x_sub_1);
  auto t3 = ADD(x, SUB(x, c1));
  auto t4 = ADD(x, x_and_1);
  auto t5 = ADD(y, x_and_1);
  auto t6 = ADD(y, BVAND(x, c1));
  auto tc1 = L_AND(L_AND(L_AND(BVULT(x, c8), BVULT(c4, x)),
                         L_AND(BVULT(y, c4), BVULT(c1, y))),
                   EQ(BVAND(x, y), c0));
  auto tc2 = L_AND(BVULT(c2, x), BVULT(c2, y));
  auto tc3 = L_AND(L_AND(BVULT(c2, x), BVULT(c2, y)), BVULT(x, y));

  auto e1 = ITE(EQ(SUB(x, y), z), ADD(a, b), c);
  auto e2 = ITE(EQ(SUB(x, y), x), ADD(a, b), SUB(a, c));
  auto e3 = ITE(EQ(BVAND(x, y), x), BVAND(a, b), BVOR(a, c));
  auto e4a = ITE(EQ(SUB(u, SUB(v, w)), d), ADD(e, f), g);
  auto e4b = ITE(
      EQ(EQ(BVAND(a, SUB(e1, e2)), f), EQ(BVAND(x, y), x)), ADD(e4a, f), e1);

  PartialModelGen pt(s);

  CheckPartialModel_bitlevel(EQ(x_plus_1, u));
  CheckPartialModel_bitlevel(EQ(x_sub_1, u));
  CheckPartialModel_bitlevel(EQ(x_and_1, u));
  CheckPartialModel_bitlevel(EQ(x_plus_y, u));
  CheckPartialModel_bitlevel(EQ(x_sub_y, u));
  CheckPartialModel_bitlevel(EQ(x_and_y, u));
  CheckPartialModel_bitlevel(EQ(t0, u));
  CheckPartialModel_bitlevel(EQ(t1, u));
  CheckPartialModel_bitlevel(EQ(t2, u));
  CheckPartialModel_bitlevel(EQ(t3, u));
  CheckPartialModel_bitlevel(EQ(t4, u));
  CheckPartialModel_bitlevel(EQ(t5, u));
  CheckPartialModel_bitlevel(EQ(t6, u));
  CheckPartialModel_bitlevel(EQ(e1, u));
  CheckPartialModel_bitlevel(EQ(e2, u));
  CheckPartialModel_bitlevel(EQ(e3, u));
  CheckPartialModel_bitlevel(EQ(e4a, u));
  CheckPartialModel_bitlevel(EQ(e4b, u));
  CheckPartialModel_bitlevel(tc1);
  CheckPartialModel_bitlevel(tc2);
  CheckPartialModel_bitlevel(tc3);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedDynamicCoiUnitTests,
                         DynamicCoiUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
