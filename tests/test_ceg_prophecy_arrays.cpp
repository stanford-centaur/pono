#include <tuple>

#include "core/prop.h"
#include "core/rts.h"
#include "engines/ceg_prophecy_arrays.h"
#include "engines/cegar_ops_uf.h"
#include "engines/ic3ia.h"
#include "gtest/gtest.h"
#include "options/options.h"
#include "smt-switch/smt.h"
#include "smt/available_solvers.h"
#include "utils/make_provers.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class CegProphecyArraysTest
    : public ::testing::TestWithParam<tuple<SolverEnum, SolverEnum>>
{
 protected:
  void SetUp() override
  {
    opts.smt_solver_ = get<0>(GetParam());
    opts.smt_interpolator_ = get<1>(GetParam());
    s = create_solver(opts.smt_solver_);
    s->set_opt("produce-unsat-assumptions", "true");
  }
  PonoOptions opts;
  SmtSolver s;
};

TEST_P(CegProphecyArraysTest, Simple)
{
  RelationalTransitionSystem rts(s);
  Sort intsort = rts.make_sort(INT);
  Sort arrsort = rts.make_sort(ARRAY, intsort, intsort);
  Term i = rts.make_statevar("i", intsort);
  Term j = rts.make_statevar("j", intsort);
  Term d = rts.make_statevar("d", intsort);
  Term a = rts.make_statevar("a", arrsort);

  Term constarr0 = rts.make_term(rts.make_term(0, intsort), arrsort);
  rts.set_init(rts.make_term(Equal, a, constarr0));
  rts.assign_next(
      a,
      rts.make_term(Ite,
                    rts.make_term(Lt, d, rts.make_term(200, intsort)),
                    rts.make_term(Store, a, i, d),
                    a));

  Term prop_term = rts.make_term(
      Lt, rts.make_term(Select, a, j), rts.make_term(200, intsort));
  SafetyProperty prop(s, prop_term);
  std::shared_ptr<SafetyProver> cegp =
      make_ceg_proph_prover(INTERP, prop, rts, s, opts);
  ProverResult r = cegp->check_until(5);
  ASSERT_EQ(r, ProverResult::TRUE);
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedCegProphecyArraysTest,
    CegProphecyArraysTest,
    testing::Combine(
        testing::ValuesIn(filter_solver_enums({ THEORY_INT, smt::CONSTARR })),
        testing::ValuesIn(filter_interpolator_enums({ THEORY_INT }))));

class CegProphecyArraysBVTest : public CegProphecyArraysTest
{
};

// the prover built for --engine ic3ia --ceg-bv-arith --ceg-prophecy-arrays,
// exposing the state of the wrapped CegProphecyArrays after initialize()
class CegarOpsUfCegpIC3IA : public CegarOpsUf<CegProphecyArrays<IC3IA>>
{
 public:
  using CegarOpsUf<CegProphecyArrays<IC3IA>>::CegarOpsUf;

  bool has_finite_index_sort() const { return has_finite_index_sort_; }
  const TransitionSystem & abs_ts() const { return abs_ts_; }
};

// CegarOpsUf fills the system of the CegProphecyArrays it wraps only when it
// initializes, so the index sorts must be checked then
TEST_P(CegProphecyArraysBVTest, FiniteIndexSortUnderCegarOpsUf)
{
  if (opts.smt_solver_ == BTOR) {
    GTEST_SKIP() << "Boolector cannot declare the uninterpreted sort that CEGP "
                    "abstracts arrays to";
  }
  RelationalTransitionSystem rts(s);
  Sort bvsort = rts.make_sort(BV, 8);
  Sort arrsort = rts.make_sort(ARRAY, bvsort, bvsort);
  Term i = rts.make_statevar("i", bvsort);
  Term a = rts.make_statevar("a", arrsort);
  rts.assign_next(i, rts.make_term(BVAdd, i, rts.make_term(1, bvsort)));
  rts.assign_next(a, rts.make_term(Store, a, i, i));
  SafetyProperty prop(s, rts.make_term(Equal, rts.make_term(Select, a, i), i));

  CegarOpsUfCegpIC3IA prover(prop, rts, s, opts);
  prover.set_ops_to_abstract({ BVAdd });
  prover.initialize();
  EXPECT_TRUE(prover.has_finite_index_sort());
  // so no lambda indices, which are integers, are created
  for (const auto & sv : prover.abs_ts().statevars()) {
    EXPECT_NE(sv->get_sort()->get_sort_kind(), INT);
  }
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedCegProphecyArraysBVTest,
    CegProphecyArraysBVTest,
    testing::Combine(
        testing::ValuesIn(filter_solver_enums({ THEORY_BV })),
        testing::ValuesIn(filter_interpolator_enums({ THEORY_BV }))));

}  // namespace pono_tests
