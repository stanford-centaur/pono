/*!
 * \file tableau.cpp
 * \brief Pure tableau/latch-building primitives for SVA/LTL encoding.
 * \author Áron Ricardo Perez-Lopez
 * \date 2026
 * \copyright See the LICENSE file in the top-level source directory.
 *
 * See tableau.h for what this class covers and why it holds no reference
 * back to SystemVerilogEncoder.
 *
 * LTL tableau
 * -----------
 * Properties that are not pure safety are translated with a standard
 * symbolic LTL tableau (temporal testers).  Every temporal operator
 * introduces a one-step "promise" tester:
 *
 *   * a free 1-bit input  n  guessing whether the operator's body
 *     holds at the *next* cycle (this is the operator's value now);
 *   * a 1-bit state  s  with  s' = n  (it remembers the guess), plus
 *     the per-cycle consistency constraint  s == body, which forces
 *     the guess made one cycle earlier to have been correct.
 *
 * Greatest-fixpoint operators (G, R) need no fairness.  Least-fixpoint
 * operators (F, strong-U) additionally emit a justice (GF) condition
 * that discharges the eventuality, ruling out lassos where a promise
 * is made forever but never fulfilled.  Callers (ltl_to_sat() in
 * assertion_walker.cpp) push negation to the leaves on the fly, so the
 * testers built here always match the negation-normal form of the
 * (negated) property.
 */
#include "frontends/systemverilog/tableau.h"

using namespace smt;
using namespace std;

namespace pono {

Tableau::Tableau(FunctionalTransitionSystem & fts,
                 const smt::SmtSolver & solver)
    : fts_(fts), solver_(solver)
{
}

string Tableau::make_name(const string & name_prefix, const string & name) const
{
  if (name_prefix.empty()) return name;
  return name_prefix + "." + name;
}

Term Tableau::make_history_chain(const Term & value,
                                 uint32_t n,
                                 const string & name_prefix)
{
  Sort sort = value->get_sort();
  Term zero = solver_->make_term(0, sort);
  Term link = value;
  for (uint32_t i = 0; i < n; ++i) {
    Term latch = fts_.make_statevar(
        make_name(name_prefix,
                  "__sva_past_" + std::to_string(latch_counter_++)),
        sort);
    fts_.constrain_init(solver_->make_term(Equal, latch, zero));
    fts_.assign_next(latch, link);
    link = latch;
  }
  return link;
}

Term Tableau::delay_bool(const Term & cond,
                         uint32_t n,
                         const string & name_prefix)
{
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  Term cond_bv = solver_->make_term(Ite, cond, one_bv1, zero_bv1);
  Term delayed_bv = make_history_chain(cond_bv, n, name_prefix);
  return solver_->make_term(Equal, delayed_bv, one_bv1);
}

smt::Term Tableau::init_flag(const string & name_prefix)
{
  if (ltl_init_flag_) return ltl_init_flag_;
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  Term flag =
      fts_.make_statevar(make_name(name_prefix, "__ltl_init_flag"), bv1);
  fts_.constrain_init(solver_->make_term(Equal, flag, one_bv1));
  fts_.assign_next(flag, zero_bv1);
  ltl_init_flag_ = solver_->make_term(Equal, flag, one_bv1);
  return ltl_init_flag_;
}

smt::Term Tableau::before_cycle(uint32_t k, const string & name_prefix)
{
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  // "is cycle 0" pulse, delayed by 1..k-1 steps gives "is cycle i" for
  // each i in [1, k); their disjunction (with the undelayed pulse) is
  // true exactly during cycles [0, k).
  Term pulse = init_flag(name_prefix);
  Term result = pulse;
  Term link = solver_->make_term(Ite, pulse, one_bv1, zero_bv1);
  for (uint32_t i = 1; i < k; ++i) {
    Term latch = fts_.make_statevar(
        make_name(name_prefix,
                  "__before_cycle_" + std::to_string(latch_counter_++)),
        bv1);
    fts_.constrain_init(solver_->make_term(Equal, latch, zero_bv1));
    fts_.assign_next(latch, link);
    result = solver_->make_term(
        Or, result, solver_->make_term(Equal, latch, one_bv1));
    link = latch;
  }
  return result;
}

smt::Term Tableau::disable_window(const Term & disable_cond,
                                  uint32_t window,
                                  const string & name_prefix)
{
  if (!disable_cond) return Term();
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  // OR the disable condition together with its 1..window-cycle-delayed
  // versions, so a `disable iff` that was true anywhere in an
  // antecedent's shift window still exempts the whole property, not
  // just the single cycle the check ends up anchored at.
  Term result = disable_cond;
  Term link = solver_->make_term(Ite, disable_cond, one_bv1, zero_bv1);
  for (uint32_t i = 0; i < window; ++i) {
    Term latch = fts_.make_statevar(
        make_name(name_prefix,
                  "__disable_hist_" + std::to_string(latch_counter_++)),
        bv1);
    fts_.constrain_init(solver_->make_term(Equal, latch, zero_bv1));
    fts_.assign_next(latch, link);
    result = solver_->make_term(
        Or, result, solver_->make_term(Equal, latch, one_bv1));
    link = latch;
  }
  return result;
}

smt::Term Tableau::make_X(const smt::Term & phi, const string & name_prefix)
{
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  uint32_t id = latch_counter_++;
  Term n = fts_.make_inputvar(
      make_name(name_prefix, "__ltl_n_" + std::to_string(id)), bv1);
  Term s = fts_.make_statevar(
      make_name(name_prefix, "__ltl_s_" + std::to_string(id)), bv1);
  fts_.assign_next(s, n);  // s@(t+1) = n@t  (remember the guess)
  // The guess made at t-1 about phi@t must have been correct.  The
  // cycle-0 instance only pins the otherwise-unused initial s.
  Term phi_bv = solver_->make_term(Ite, phi, one_bv1, zero_bv1);
  fts_.add_constraint(solver_->make_term(Equal, s, phi_bv),
                      /*to_init_and_next=*/false);
  return solver_->make_term(Equal, n, one_bv1);  // sat(X phi) = guess
}

smt::Term Tableau::make_G(const smt::Term & phi, const string & name_prefix)
{
  // G phi == phi && X(G phi)
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  uint32_t id = latch_counter_++;
  Term n = fts_.make_inputvar(
      make_name(name_prefix, "__ltl_n_" + std::to_string(id)), bv1);
  Term s = fts_.make_statevar(
      make_name(name_prefix, "__ltl_s_" + std::to_string(id)), bv1);
  fts_.assign_next(s, n);
  Term n_bool = solver_->make_term(Equal, n, one_bv1);
  Term body = solver_->make_term(And, phi, n_bool);  // sat(G phi)
  Term body_bv = solver_->make_term(Ite, body, one_bv1, zero_bv1);
  fts_.add_constraint(solver_->make_term(Equal, s, body_bv),
                      /*to_init_and_next=*/false);
  return body;
}

smt::Term Tableau::make_F(const smt::Term & phi,
                          smt::TermVec & justice,
                          const string & name_prefix)
{
  // F phi == phi || X(F phi)
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  uint32_t id = latch_counter_++;
  Term n = fts_.make_inputvar(
      make_name(name_prefix, "__ltl_n_" + std::to_string(id)), bv1);
  Term s = fts_.make_statevar(
      make_name(name_prefix, "__ltl_s_" + std::to_string(id)), bv1);
  fts_.assign_next(s, n);
  Term n_bool = solver_->make_term(Equal, n, one_bv1);
  Term body = solver_->make_term(Or, phi, n_bool);  // sat(F phi)
  Term body_bv = solver_->make_term(Ite, body, one_bv1, zero_bv1);
  fts_.add_constraint(solver_->make_term(Equal, s, body_bv),
                      /*to_init_and_next=*/false);
  // Discharge: infinitely often phi holds or we stop promising F.
  justice.push_back(
      solver_->make_term(Or, phi, solver_->make_term(Not, n_bool)));
  return body;
}

smt::Term Tableau::make_R(const smt::Term & a,
                          const smt::Term & b,
                          const string & name_prefix)
{
  // a R b == b && (a || X(a R b))   (greatest fixpoint, no fairness)
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  uint32_t id = latch_counter_++;
  Term n = fts_.make_inputvar(
      make_name(name_prefix, "__ltl_n_" + std::to_string(id)), bv1);
  Term s = fts_.make_statevar(
      make_name(name_prefix, "__ltl_s_" + std::to_string(id)), bv1);
  fts_.assign_next(s, n);
  Term n_bool = solver_->make_term(Equal, n, one_bv1);
  Term body = solver_->make_term(
      And, b, solver_->make_term(Or, a, n_bool));  // sat(a R b)
  Term body_bv = solver_->make_term(Ite, body, one_bv1, zero_bv1);
  fts_.add_constraint(solver_->make_term(Equal, s, body_bv),
                      /*to_init_and_next=*/false);
  return body;
}

smt::Term Tableau::make_U(const smt::Term & a,
                          const smt::Term & b,
                          smt::TermVec & justice,
                          const string & name_prefix)
{
  // a U b (strong) == b || (a && X(a U b))
  Sort bv1 = solver_->make_sort(BV, 1);
  Term one_bv1 = solver_->make_term(1, bv1);
  Term zero_bv1 = solver_->make_term(0, bv1);
  uint32_t id = latch_counter_++;
  Term n = fts_.make_inputvar(
      make_name(name_prefix, "__ltl_n_" + std::to_string(id)), bv1);
  Term s = fts_.make_statevar(
      make_name(name_prefix, "__ltl_s_" + std::to_string(id)), bv1);
  fts_.assign_next(s, n);
  Term n_bool = solver_->make_term(Equal, n, one_bv1);
  Term body = solver_->make_term(
      Or, b, solver_->make_term(And, a, n_bool));  // sat(a U b)
  Term body_bv = solver_->make_term(Ite, body, one_bv1, zero_bv1);
  fts_.add_constraint(solver_->make_term(Equal, s, body_bv),
                      /*to_init_and_next=*/false);
  // Discharge: infinitely often b holds or a U b is false.
  justice.push_back(solver_->make_term(Or, solver_->make_term(Not, body), b));
  return body;
}

}  // namespace pono
