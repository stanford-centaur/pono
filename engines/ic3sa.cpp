/*********************                                                  */
/*! \file ic3sa.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief IC3 with Syntax-Guided Abstraction based on
**
**        Model Checking of Verilog RTL Using IC3 with Syntax-Guided
*Abstraction.
**            -- Aman Goel, Karem Sakallah
**
**
**  within Pono, we are building on the bit-level IC3 instead of directly
**  on IC3Base, because a lot of the functionality is the same
**  In particular, we don't need to override inductive generalization
**
**/

#include "engines/ic3sa.h"

#include "assert.h"
#include "utils/term_walkers.h"

using namespace smt;
using namespace std;

namespace pono {

// main IC3SA implementation

IC3SA::IC3SA(Property & p, smt::SolverEnum se) : super(p, se) {}

IC3SA::IC3SA(Property & p, const smt::SmtSolver & s) : super(p, s) {}

IC3SA::IC3SA(const PonoOptions & opt, Property & p, smt::SolverEnum se)
    : super(opt, p, se)
{
}

IC3SA::IC3SA(const PonoOptions & opt, Property & p, const smt::SmtSolver & s)
    : super(opt, p, s)
{
}

IC3Formula IC3SA::get_model_ic3formula(TermVec * out_inputs,
                                       TermVec * out_nexts) const
{
  EquivalenceClasses ec = get_equivalence_classes_from_model();

  // now create a cube expressing this partition
  TermVec cube_lits;
  for (const auto & sortelem : ec) {
    const Sort & sort = sortelem.first;

    // TODO: play around with heuristics for the representative
    //       to add disequalities over
    //       e.g. we're not adding all possible disequalities,
    //       just choosing a representative from each equivalence
    //       class and adding a disequality to encode the distinctness

    //       currently preferring symbol > generic term > value

    // representatives of the different classes of this sort
    TermVec representatives;
    for (const auto & elem : sortelem.second) {
      const Term & val = elem.first;
      assert(val->is_value());

      const UnorderedTermSet & terms = elem.second;
      Term repr = val;
      bool found_repr = false;
      bool repr_val = true;

      UnorderedTermSet::const_iterator end = terms.cend();
      UnorderedTermSet::const_iterator it = terms.cbegin();
      Term last = *it;

      while (it != end) {
        const Term & term = *(++it);
        assert(last->get_sort() == term->get_sort());
        cube_lits.push_back(solver_->make_term(Equal, last, term));
        last = term;

        // TODO: see if a DisjointSet would make this easier
        // update representative for this class
        if (!found_repr) {
          if (term->is_symbolic_const()) {
            repr = term;
            repr_val = false;
            found_repr = true;
          } else if (!term->is_value() && repr_val) {
            repr = term;
            repr_val = false;
          }
        }
      }
    }

    // add disequalities between each pair of representatives from
    // different equivalent classes
    for (size_t i = 0; i < representatives.size(); ++i) {
      for (size_t j = i + 1; j < representatives.size(); ++j) {
        const Term & ti = representatives.at(i);
        const Term & tj = representatives.at(j);
        // should never get the same representative term from different classes
        assert(ti != tj);
        cube_lits.push_back(solver_->make_term(Distinct, ti, tj));
      }
    }
  }

  IC3Formula cube = ic3formula_conjunction(cube_lits);
  assert(ic3formula_check_valid(cube));
  return cube;
}

bool IC3SA::ic3formula_check_valid(const IC3Formula & u) const
{
  throw PonoException("IC3SA::ic3formula_check_valid NYI");
}

IC3Formula IC3SA::generalize_predecessor(size_t i, const IC3Formula & c)
{
  throw PonoException("IC3SA::generalize_predecessor NYI");
}

void IC3SA::check_ts() const
{
  // TODO: add support for arrays

  for (const auto & sv : ts_->statevars())
  {
    SortKind sk = sv->get_sort()->get_sort_kind();
    if (sk != BOOL && sk != BV)
    {
      throw PonoException("IC3SA currently only supports bit-vectors");
    }
  }
  for (const auto & iv : ts_->inputvars())
  {
    SortKind sk = iv->get_sort()->get_sort_kind();
    if (sk != BOOL && sk != BV)
    {
      throw PonoException("IC3SA currently only supports bit-vectors");
    }
  }
}

RefineResult IC3SA::refine() { throw PonoException("IC3SA::refine NYI"); }

void IC3SA::initialize()
{
  super::initialize();

  // set up initial term abstraction by getting all subterms
  // TODO consider starting with only a subset -- e.g. variables
  SubTermCollector stc(solver_, false);
  stc.collect_subterms(ts_->init());
  stc.collect_subterms(ts_->trans());
  stc.collect_subterms(bad_);
  term_abstraction_ = stc.get_subterms();

  throw PonoException("IC3SA::initialize not completed");
}

// IC3SA specific methods

EquivalenceClasses IC3SA::get_equivalence_classes_from_model() const
{
  // assumes the solver state is sat
  EquivalenceClasses ec;
  for (auto elem : term_abstraction_) {
    const Sort & sort = elem.first;
    const UnorderedTermSet & terms = elem.second;

    // TODO figure out if a DisjointSet is a better data structure
    //      will need to keep track of all terms in each partition though

    ec[sort] = std::unordered_map<smt::Term, smt::UnorderedTermSet>();
    std::unordered_map<smt::Term, smt::UnorderedTermSet> & m = ec.at(sort);
    for (auto t : terms) {
      Term val = solver_->get_value(t);
      m[val].insert(t);
    }
  }
  return ec;
}

}  // namespace pono
