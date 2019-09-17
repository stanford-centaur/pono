#pragma once

#include "rts.h"
#include "prop.h"
#include "unroller.h"

#include "smt-switch/smt.h"
#include "proverresult.h"

namespace cosa
{

class KInduction
{
public:
  KInduction(const Property &p, smt::SmtSolver &solver);
  ~KInduction();

  void initialize();

  ProverResult check_until(int k);

  ProverResult prove();

  bool witness(std::vector<smt::UnorderedTermMap> &out);

private:

  bool base_step(int i);
  bool inductive_step(int i);

  void add_simple_path_constraint(int i, int j);

  const RelationalTransitionSystem &ts_;
  const Property &property_;

  smt::SmtSolver &solver_;
  Unroller unroller_;

  int reached_k_;

  smt::Term simple_path_;

}; // class KInduction

} // namespace cosa
