#ifndef PROP_H
#define PROP_H

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

private:

  const RelationalTransitionSystem &ts_;
  
  smt::Term prop_;
  
}; // class Property

} // namespace cosa

#endif
