#pragma once

#include <stdio.h>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>

#include "assert.h"
#include "exceptions.h"
#include "smt-switch/smt.h"

/* node to help construct syntax tree in the parser */
namespace cosa {
struct SMVnode
{
  enum Nodetype
  {
    TermNode,
    SortNode
  };
  smt::Term tm;
  smt::Sort st;
  Nodetype nt;
  SMVnode(smt::Sort s)
  {
    nt = SortNode;
    st = s;
  }
  SMVnode(smt::Term t)
  {
    nt = TermNode;
    tm = t;
  }
  Nodetype getNodeType() { return nt; }
  smt::Sort getSort() { return st; }
  smt::Term getTerm() { return tm; }
};  // struct SMVNode
}  // namespace cosa