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
  enum BVtype{
    Signed,
    Unsigned,
    BVnot
  };
  smt::Term tm;
  smt::Sort st;
  Nodetype nt;
  BVtype bvt;
  SMVnode(smt::Sort s)
  {
    nt = SortNode;
    st = s;
    bvt = BVnot;
  }
  SMVnode(smt::Term t)
  {
    nt = TermNode;
    tm = t;
    bvt = BVnot;
  }
  SMVnode(smt::Sort s, BVtype bvtype)
  {
    nt = SortNode;
    st = s;
    bvt = bvtype;
  }
  SMVnode(smt::Term t, BVtype bvtype)
  {
    nt = TermNode;
    tm = t;
    bvt = bvtype;
  }
  Nodetype getNodeType() { return nt; }
  BVtype getBVType() { return bvt; }
  smt::Sort getSort() { return st; }
  smt::Term getTerm() { return tm; }
};  // struct SMVNode
}  // namespace cosa