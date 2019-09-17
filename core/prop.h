#pragma once

#include "smt-switch/smt.h"
#include "rts.h"

namespace cosa
{

class Property
{
public:
  Property(const RelationalTransitionSystem &ts, smt::Term p);
  ~Property();

  const smt::Term prop() const { return prop_; }

  const RelationalTransitionSystem &transition_system() const {
    return ts_;
  }

private:

  const RelationalTransitionSystem &ts_;
  
  smt::Term prop_;
  
}; // class Property

} // namespace cosa
