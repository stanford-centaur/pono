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
  enum Type
  {
    Signed,
    Unsigned,
    Integer,
    Boolean,
    BV,
    BVnot
  };
  smt::Term tm;
  smt::Sort st;
  Nodetype nt;
  Type bvt;
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
  SMVnode(smt::Sort s, Type type)
  {
    nt = SortNode;
    st = s;
    bvt = type;
  }
  SMVnode(smt::Term t, Type type)
  {
    nt = TermNode;
    tm = t;
    bvt = type;
  }
  Nodetype getNodeType() { return nt; }
  Type getType() { return bvt; }
  smt::Sort getSort() { return st; }
  smt::Term getTerm() { return tm; }
};  // struct SMVNode
}  // namespace cosa