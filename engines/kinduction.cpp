#include "kinduction.h"
#include "utils/logger.h"

using namespace smt;

namespace cosa {

KInduction::KInduction(const Property & p, smt::SmtSolver & solver)
    : ts_(p.transition_system()),
      property_(p),
      solver_(solver),
      unroller_(ts_, solver_)
{
  initialize();
}

KInduction::~KInduction() {}

void KInduction::initialize()
{
  reached_k_ = -1;
  // NOTE: There's an implicit assumption that this solver is only used for
  // model checking once Otherwise there could be conflicting assertions to
  // the solver or it could just be polluted with redundant assertions in the
  // future we can use solver_->reset_assertions(), but it is not currently
  // supported in boolector
  simple_path_ = solver_->make_value(true);
}

ProverResult KInduction::check_until(int k)
{
  for (int i = 0; i <= k; ++i)
  {
    logger.log(1, "Checking k-induction base case at bound: {}", i);
    if (!base_step(i))
    {
      return ProverResult::FALSE;
    }
    logger.log(1, "Checking k-induction inductive step at bound: {}", i);
    if (inductive_step(i))
    {
      return ProverResult::TRUE;
    }
  }
  return ProverResult::UNKNOWN;
}

ProverResult KInduction::prove()
{
  for (int i = 0;; ++i)
  {
    if (!base_step(i))
    {
      return ProverResult::FALSE;
    }
    if (inductive_step(i))
    {
      return ProverResult::TRUE;
    }
  }
  return ProverResult::UNKNOWN;
}

bool KInduction::witness(std::vector<UnorderedTermMap> & out)
{
  // TODO: make sure the solver state is SAT

  for (int i = 0; i <= reached_k_; ++i)
  {
    out.push_back(UnorderedTermMap());
    UnorderedTermMap & map = out.back();

    for (auto v : ts_.states())
    {
      Term vi = unroller_.at_time(v, i);
      Term r = solver_->get_value(vi);
      map[v] = r;
    }

    for (auto v : ts_.inputs())
    {
      Term vi = unroller_.at_time(v, i);
      Term r = solver_->get_value(vi);
      map[v] = r;
    }
  }

  return true;
}

bool KInduction::base_step(int i)
{
  if (i <= reached_k_)
  {
    return true;
  }

  Term bad = solver_->make_term(PrimOp::Not, property_.prop());

  solver_->push();
  solver_->assert_formula(unroller_.at_time(ts_.init(), 0));
  solver_->assert_formula(unroller_.at_time(bad, i));
  Result r = solver_->check_sat();
  if (r.is_sat())
  {
    ++reached_k_;
    return false;
  }
  solver_->pop();

  solver_->assert_formula(unroller_.at_time(ts_.trans(), i));
  solver_->assert_formula(unroller_.at_time(property_.prop(), i));

  return true;
}

bool KInduction::inductive_step(int i)
{
  if (i <= reached_k_)
  {
    return false;
  }

  Term bad = solver_->make_term(PrimOp::Not, property_.prop());

  solver_->push();
  solver_->assert_formula(simple_path_);
  solver_->assert_formula(unroller_.at_time(bad, i + 1));
  Result r = solver_->check_sat();

  Term f = solver_->make_value(false);
  Term constraint;
  bool keep_checking = true;
  bool added_to_simple_path = false;

  while (keep_checking)
  {
    if (r.is_unsat())
    {
      return true;
    }

    added_to_simple_path = false;
    for (int j = 0; j < i && !added_to_simple_path; ++j)
    {
      for (int l = j + 1; l < i; ++l)
      {
        constraint = simple_path_constraint(j, l);
        if (solver_->get_value(constraint) == f)
        {
          simple_path_ =
              solver_->make_term(PrimOp::And, simple_path_, constraint);
          added_to_simple_path = true;
          break;
        }
      }
    }

    if (added_to_simple_path)
    {
      solver_->assert_formula(constraint);
      r = solver_->check_sat();
      keep_checking = (r.is_sat() && added_to_simple_path);
    }
    else
    {
      keep_checking = false;
    }
  }

  solver_->pop();

  ++reached_k_;

  return false;
}

Term KInduction::simple_path_constraint(int i, int j)
{
  // TODO: what if there are no states?
  //       kind of a weird situation, but possible -- don't want to assume false
  Term disj = solver_->make_value(false);
  for (auto v : ts_.states())
  {
    Term vi = unroller_.at_time(v, i);
    Term vj = unroller_.at_time(v, j);
    Term eq = solver_->make_term(PrimOp::Equal, vi, vj);
    Term neq = solver_->make_term(PrimOp::Not, eq);
    disj = solver_->make_term(PrimOp::Or, disj, neq);
  }
  return disj;
}

}  // namespace cosa
