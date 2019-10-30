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

 protected:
  bool base_step(int i);
  bool inductive_step(int i);

  smt::Term simple_path_constraint(int i, int j);
  bool check_simple_path_lazy(int i);

  const RelationalTransitionSystem & ts_;
  const Property & property_;

  smt::SmtSolver & solver_;
  Unroller unroller_;

  int reached_k_;

  smt::Term init0_;
  smt::Term bad_;
  smt::Term false_;
  smt::Term simple_path_;

};  // class KInduction

}  // namespace cosa
