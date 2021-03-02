/*********************                                                  */
/*! \file partial_model.h
** \verbatim
** Top contributors (to current version):
**   Hongce Zhang
** This file is part of the pono project.
** Copyright (c) 2020 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Class for performing dynamic cone of influence reduction based
**        on the model from the solver. This is essentially extracting a
**        partial model. (see tests/test_partial_model.cpp for usage)
**
**/

#pragma once

#include "utils/sygus_ic3formula_helper.h"
#include "engines/ic3base.h"

namespace pono {

class PartialModelGen {
public:
  /** This class computes the cone of influence on construction
   *  The current implementation does not have internal cache,
   *  but in the future maybe we can cache some of the results
   *  @param the solver where the assertions were made
   */
  PartialModelGen(smt::SmtSolver & solver) : solver_(solver) { }
    
  // disallow copy construct/assign
  PartialModelGen(const PartialModelGen &) = delete;
  PartialModelGen & operator=(const PartialModelGen &) = delete;
  
protected:
  // let's keep a reference to the solver since we need to add terms
  smt::SmtSolver & solver_;

  // for the DFS, will not use the stack but use one reference here
  std::unordered_set<smt::Term> dfs_walked_;
  std::unordered_set<smt::Term> dfs_vars_;
  void dfs_walk(const smt::Term & ast);

  // conditon var buffer
  void GetVarList(const smt::Term & ast);

public:

  /** This class computes the variables that need to
   *  appear in the partial model of ast
   *  @param the ast to walk
   *  @param (output) the set of variables
   */
  void GetVarList(const smt::Term & ast, 
    std::unordered_set<smt::Term> & out_vars);

  /** This class computes the variables that need to
   *  appear in the partial model of asts in the vector
   *  @param the vector of ast to walk
   *  @param (output) the set of variables
   */
  void GetVarListForAsts(const smt::TermVec & asts, 
    smt::UnorderedTermSet & out_vars);

  /** This class computes the variables that need to
   *  appear in the partial model of asts in the vector
   *  @param the ast to walk
   *  @return the partial model in ic3formula
   */
  IC3Formula GetPartialModel(const smt::Term & ast);

  /** This class computes the variables that need to
   *  appear in the partial model of asts in the vector
   *  @param the ast to walk
   *  @return the partial model and the var/val cube
   */
  std::pair<IC3Formula,syntax_analysis::IC3FormulaModel> 
    GetPartialModelInCube(const smt::Term & ast);


  // add an API to use buffers 
};

}  // namespace pono
