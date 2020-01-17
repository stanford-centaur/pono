#include "fts.h"
#include "smt-switch/smt.h"

// based on implementation here: https://github.com/cristian-mattarei/CoSA/blob/master/cosa/utils/formula_mngm.py

namespace cosa
{
  // symbolically sorts a list of boolean terms
  // true before false
  class BoolSortingNetwork
  {
  public:
    BoolSortingNetwork(smt::SmtSolver solver, smt::TermVec boolvec) : solver_(solver), boolvec_(boolvec)
    {
      do_sorting();
    }

    smt::TermVec sorted_boolvec() { return sorted_boolvec_; };

  protected:
    void do_sorting();
    smt::TermVec sorting_network_helper(smt::TermVec boolvec) const;
    smt::TermVec two_comparator(smt::Term bool0, smt::Term bool1) const;
    smt::TermVec merge(smt::TermVec boolvec0, smt::TermVec boolvec1) const;

    smt::SmtSolver solver_;
    smt::TermVec sorted_boolvec_;
    smt::TermVec boolvec_;
  };
}
