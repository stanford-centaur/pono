#pragma once

#include "prop.h"
#include "prover.h"
#include "rts.h"
#include "unroller.h"

#include "proverresult.h"
#include "smt-switch/smt.h"

namespace cosa {

class KInduction : public Prover
{
 public:
  KInduction(const Property & p, smt::SmtSolver & solver);
  ~KInduction();

  void initialize();

  ProverResult check_until(int k);

  ProverResult prove();

  bool witness(std::vector<smt::UnorderedTermMap> & out);

 private:
  bool base_step(int i);
  bool inductive_step(int i);

  void add_simple_path_constraint(int i, int j);

  const RelationalTransitionSystem & ts_;
  const Property & property_;

  smt::SmtSolver & solver_;
  Unroller unroller_;

  int reached_k_;

  smt::Term simple_path_;

};  // class KInduction

}  // namespace cosa
