#pragma once

#include <stdio.h>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include "assert.h"
#include "utils/exceptions.h"
#include "smt-switch/smt.h"

/* node to help construct syntax tree in the parser */
namespace pono {
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
    Real,
    Boolean,
    WordArray,
    BV,
    Default
  };
  smt::Term tm;
  smt::Sort st;
  Nodetype nt;
  Type bvt;
  SMVnode(smt::Sort s)
  {
    nt = SortNode;
    st = s;
    bvt = Default;
  }
  SMVnode(smt::Term t)
  {
    nt = TermNode;
    tm = t;
    bvt = Default;
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
}  // namespace pono
