#include "prop.h"
#include "exceptions.h"

using namespace smt;

namespace cosa
{

  Properties::Properties(const RelationalTransitionSystem &ts):
    ts_(ts)
  {
  }

  Properties::~Properties()
  {
  }

  void Properties::add_prop(const Term p)
  {
    prop_.push_back(p);
  }

  Term Properties::get_prop(unsigned int i)
  {
    if (i >= prop_.size()) {
      throw CosaException("Property index is out of bound");
    }
    
    return prop_[i];
  }
  
  const std::vector<Term> &Properties::prop()
  {
    return prop_;
  }
  
} // namespace cosa
