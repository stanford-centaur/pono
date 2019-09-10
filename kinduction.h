#pragma once

#include "rts.h"
#include "prop.h"
#include "unroller.h"

#include "smt-switch/smt.h"

namespace cosa
{

class KInduction
{
public:
  KInduction(const Property &p, smt::SmtSolver &solver);
  ~KInduction();

  void initialize();

  bool check_until(size_t k);
  
private:

  bool base_step(size_t i);
  bool inductive_step(size_t i);

  void add_simple_path_constraint(size_t i, size_t j);

  const RelationalTransitionSystem &ts_;
  const Property &property_;

  smt::SmtSolver &solver_;
  Unroller unroller_;

  size_t reached_k_;
  
}; // class KInduction
  
} // namespace cosa
