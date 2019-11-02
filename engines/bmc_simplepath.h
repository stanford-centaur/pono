#pragma once

#include "kinduction.h"

namespace cosa {

class BmcSimplePath : public KInduction
{
 public:
  BmcSimplePath(const Property & p, smt::SmtSolver & solver);
  ~BmcSimplePath();

  typedef KInduction super;

  ProverResult check_until(int k) override;
  ProverResult prove() override;

 private:
  bool cover_step(int i);

};  // BmcSimplePath

}  // namespace cosa
