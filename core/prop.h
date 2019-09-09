#ifndef PROP_H
#define PROP_H

#include "smt-switch/smt.h"
#include "rts.h"

namespace cosa
{

class Properties
{
public:
  Properties(const RelationalTransitionSystem &ts);
  ~Properties();

  void add_prop(const smt::Term p);

  smt::Term get_prop(unsigned int i);
  
  const std::vector<smt::Term> &prop();

private:

  const RelationalTransitionSystem &ts_;
  
  std::vector<smt::Term> prop_;
  
}; // class Properties

} // namespace cosa

#endif
