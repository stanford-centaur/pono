#pragma once

#include <stdio.h>

#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include "assert.h"
#include "utils/exceptions.h"
#include "smt-switch/smt.h"

using namespace std;
/* node to help construct syntax tree in the parser */
/* SMV node is the basic node type in ast
 *  To flatten modular SMV: iterate twice before building system
 *  Two types of state variable (Vartype): basic type and module type are
 * assigned when element_node stores a vector of SMVNode * that contain the
 * corresponding type NodeMType NodeMType includes VAR, IVAR, FROZENVAR, FUN
 * DEFINE, INIT, INVAR, TRANS, ASSIGN
 */

namespace pono {
class module_node;
struct SMVnode
{
 public:
  enum Nodetype
  {
    TermNode,
    SortNode
  };
  enum Vartype
  {
    BasicT,
    ModuleT
  };
  enum NodeMtype
  {
    DEFINE,
    CONSTANTS,
    ASSIGN,
    IVAR,
    VAR,
    FROZENVAR,
    FUN,
    INIT,
    TRANS,
    INVAR,
    INVARSPEC
  };
  enum Type
  {
    Signed,
    Unsigned,
    Integer,
    Real,
    Boolean,
    // TODO shouldn't have separate WordArray and IntArray types
    // should just be Array ideally
    WordArray,
    IntArray,
    Function,
    Default
  };
  smt::Term tm;
  smt::Sort st;
  Nodetype nt;
  Type bvt;
  NodeMtype mt;
  Vartype vt;
  Type ele_type;
  SMVnode()
  {
    nt = SortNode;
    bvt = Default;
    ele_type = Default;
  }
  SMVnode(smt::Sort s)
  {
    nt = SortNode;
    st = s;
    bvt = Default;
    ele_type = Default;
  }
  SMVnode(smt::Term t)
  {
    nt = TermNode;
    tm = t;
    bvt = Default;
    ele_type = Default;
  }
  SMVnode(smt::Sort s, Type type)
  {
    nt = SortNode;
    st = s;
    bvt = type;
    ele_type = Default;
  }
  SMVnode(smt::Term t, Type type)
  {
    nt = TermNode;
    tm = t;
    bvt = type;
    ele_type = Default;
  }
  SMVnode(smt::Term t, Type type, Type element)
  {
    nt = TermNode;
    tm = t;
    bvt = type;
    ele_type = element;
  }

  Nodetype getNodeType() { return nt; }
  NodeMtype getNodeMType() { return mt; }
  Vartype getVarType() { return vt; }
  Type getType() { return bvt; }
  Type getElementType() { return ele_type; }
  smt::Sort getSort() { return st; }
  smt::Term getTerm() { return tm; }
  virtual std::string getName() { return ""; }
  virtual void generate_ostream(
      std::string name,
      std::string prefix,
      std::unordered_map<string, module_node *> module_list,
      std::unordered_map<string, string> new_prefix,
      ostream & s){};
};  // struct SMVNode

class element_node : public SMVnode
{
  std::vector<SMVnode *> pa_li;

 public:
  element_node() {}
  virtual void preprocess(
      std::string module,
      std::string prefix,
      std::unordered_map<std::string, module_node *> module_list,
      std::unordered_map<string, string> * new_prefix,
      ostream & s){};
  virtual std::vector<SMVnode *> get_list() { return pa_li; }
};

class module_node
{
  std::string module_name;
  std::string par_name;
  std::vector<string> par_li;
  unordered_map<string, SMVnode *> new_par;
  element_node * var_li;
  element_node * ivar_li;
  element_node * assign_li;
  element_node * define_li;
  element_node * frozenvar_li;
  element_node * fun_li;
  element_node * init_li;
  element_node * trans_li;
  element_node * invar_li;
  element_node * invarspec_li;

 public:
  module_node(std::string name) { module_name = name; }
  module_node(std::string name, std::vector<string> id_list)
  {
    module_name = name;
    par_li = id_list;
  }
  module_node(std::string name,
              std::vector<string> id_list,
              std::unordered_map<SMVnode::NodeMtype, element_node *> decl_map)
  {
    module_name = name;
    par_li = id_list;
    var_li = new element_node();
    ivar_li = new element_node();
    define_li = new element_node();
    assign_li = new element_node();
    frozenvar_li = new element_node();
    fun_li = new element_node();
    init_li = new element_node();
    trans_li = new element_node();
    invar_li = new element_node();
    invarspec_li = new element_node();
    for (std::unordered_map<SMVnode::NodeMtype, element_node *>::iterator it =
             decl_map.begin();
         it != decl_map.end();
         ++it) {
      if (it->first == SMVnode::DEFINE) define_li = it->second;
      if (it->first == SMVnode::ASSIGN) assign_li = it->second;
      if (it->first == SMVnode::IVAR) ivar_li = it->second;
      if (it->first == SMVnode::VAR) var_li = it->second;
      if (it->first == SMVnode::FROZENVAR) frozenvar_li = it->second;
      if (it->first == SMVnode::FUN) fun_li = it->second;
      if (it->first == SMVnode::INIT) init_li = it->second;
      if (it->first == SMVnode::TRANS) trans_li = it->second;
      if (it->first == SMVnode::INVAR) invar_li = it->second;
      if (it->first == SMVnode::INVARSPEC) invarspec_li = it->second;
    }
  }

  module_node(std::string name,
              std::unordered_map<SMVnode::NodeMtype, element_node *> decl_map)
  {
    module_name = name;
    var_li = new element_node();
    ivar_li = new element_node();
    define_li = new element_node();
    assign_li = new element_node();
    frozenvar_li = new element_node();
    fun_li = new element_node();
    init_li = new element_node();
    trans_li = new element_node();
    invar_li = new element_node();
    invarspec_li = new element_node();
    for (std::unordered_map<SMVnode::NodeMtype, element_node *>::iterator it =
             decl_map.begin();
         it != decl_map.end();
         ++it) {
      if (it->first == SMVnode::DEFINE) define_li = it->second;
      if (it->first == SMVnode::ASSIGN) assign_li = it->second;
      if (it->first == SMVnode::IVAR) ivar_li = it->second;
      if (it->first == SMVnode::VAR) var_li = it->second;
      if (it->first == SMVnode::FROZENVAR) frozenvar_li = it->second;
      if (it->first == SMVnode::FUN) fun_li = it->second;
      if (it->first == SMVnode::INIT) init_li = it->second;
      if (it->first == SMVnode::TRANS) trans_li = it->second;
      if (it->first == SMVnode::INVAR) invar_li = it->second;
      if (it->first == SMVnode::INVARSPEC) invarspec_li = it->second;
    }
  }
  unordered_map<string, SMVnode *> get_namelist() { return new_par; }
  std::string get_par() { return par_name; }
  /* process modular smv starting from main module */
  void process_main(std::unordered_map<std::string, module_node *> module_list,
           std::ostream & s);
  /* preprocess method: output to stringstream following falttened smv file format */ 
  void preprocess(std::string parent,
                  std::string prefix,
                  std::unordered_map<string, string> * new_prefix,
                  std::unordered_map<std::string, module_node *> module_list,
                  std::vector<SMVnode *> id_li,
                  std::ostream & s);
};

class type_node : public SMVnode
{
  std::string type_name;
  std::vector<SMVnode *> ex_li;
  Type ele_type;
 public:
  type_node() {}
  type_node(std::string n) { type_name = n; }
  type_node(std::string n, std::vector<SMVnode *> li)
  {
    type_name = n;
    ex_li = li;
  }
  type_node(smt::Sort s, Type type)
  {
    nt = SortNode;
    st = s;
    bvt = type;
  }

  type_node(smt::Sort s, Type arraytype, Type type)
  {
    nt = SortNode;
    st = s;
    bvt = arraytype;
    ele_type = type;
  }

  type_node(smt::Sort s)
  {
    nt = SortNode;
    st = s;
    bvt = Default;
  }
  std::string getName() { return type_name; }
  std::vector<SMVnode *> get_list() { return ex_li; }
  Type getType() { return bvt; }
  Type getElementType() { return ele_type; }
  smt::Sort getSort() { return st; }
  smt::Term getTerm() { return tm; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};

class var_node_c : public SMVnode
{
  std::string id;
  std::string type;
  pono::type_node * ty;

 public:
  var_node_c(std::string n, pono::type_node * li, Vartype vr)
  {
    ty = li;
    type = li->getName();
    id = n;
    vt = vr;
  }
  std::string getName() { return id; }
  pono::type_node * getmodtype() { return ty; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};

class ivar_node_c : public SMVnode
{
  std::string id;
  std::string type;

 public:
  ivar_node_c(std::string n, SMVnode * t)
  {
    type = t->getName();
    id = n;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};

class frozenvar_node_c : public SMVnode
{
  std::string id;
  std::string type;

 public:
  frozenvar_node_c(std::string n, SMVnode * t)
  {
    type = t->getName();
    id = n;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};

class fun_node_c : public SMVnode
{
  std::string id;
  std::string type;

 public:
  fun_node_c(std::string n, SMVnode * t)
  {
    type = t->getName();
    id = n;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};

class define_node_c : public SMVnode
{
  std::string id;
  SMVnode * ex;

 public:
  define_node_c(std::string n, SMVnode * t)
  {
    ex = t;
    id = n;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};

class assign_node_c : public SMVnode
{
  std::string pre;
  std::string id;
  SMVnode * ex;

 public:
  assign_node_c(std::string p, std::string n, SMVnode * t)
  {
    ex = t;
    id = n;
    pre = p;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};

class init_node_c : public SMVnode
{
  SMVnode * ex;

 public:
  init_node_c(SMVnode * t) { ex = t; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};

class trans_node_c : public SMVnode
{
  SMVnode * ex;

 public:
  trans_node_c(SMVnode * t) { ex = t; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};

class invar_node_c : public SMVnode
{
  SMVnode * ex;

 public:
  invar_node_c(SMVnode * t) { ex = t; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class invarspec_node_c : public SMVnode
{
  SMVnode * ex;

 public:
  invarspec_node_c(SMVnode * t) { ex = t; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};

class var_node : public element_node
{
  std::vector<var_node_c *> ex_li;

 public:
  var_node(std::vector<var_node_c *> li, NodeMtype t)
  {
    mt = t;
    ex_li = li;
  }
  void preprocess(std::string module,
                  std::string prefix,
                  std::unordered_map<std::string, module_node *> module_list,
                  std::unordered_map<string, string> * new_prefix,
                  ostream & s);

  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};

class ivar_node : public element_node
{
  std::vector<SMVnode *> ex_li;

 public:
  ivar_node(std::vector<SMVnode *> li, NodeMtype t) : ex_li(li) { mt = t; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);

  std::vector<SMVnode *> get_list() { return ex_li; }
};

class frozenvar_node : public element_node
{
  std::vector<SMVnode *> ex_li;

 public:
  frozenvar_node(std::vector<SMVnode *> li, NodeMtype t)
  {
    mt = t;
    ex_li = li;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);

  std::vector<SMVnode *> get_list() { return ex_li; }
};

class fun_node : public element_node
{
  std::vector<SMVnode *> ex_li;

 public:
  fun_node(std::vector<SMVnode *> li, NodeMtype t)
  {
    mt = t;
    ex_li = li;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);

  std::vector<SMVnode *> get_list() { return ex_li; }
};

class define_node : public element_node
{
  std::vector<SMVnode *> ex_li;

 public:
  define_node(std::vector<SMVnode *> li, NodeMtype t)
  {
    mt = t;
    ex_li = li;
  }

  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);

  std::vector<SMVnode *> get_list() { return ex_li; }
};

class assign_node : public element_node
{
  std::vector<SMVnode *> ex_li;

 public:
  assign_node(std::vector<SMVnode *> li, NodeMtype t)
  {
    mt = t;
    ex_li = li;
  }

  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);

  std::vector<SMVnode *> get_list() { return ex_li; }
};

class init_node : public element_node
{
  std::vector<SMVnode *> ex_li;

 public:
  init_node(std::vector<SMVnode *> li, NodeMtype t)
  {
    mt = t;
    ex_li = li;
  }

  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);

  std::vector<SMVnode *> get_list() { return ex_li; }
};

class trans_node : public element_node
{
  std::vector<SMVnode *> ex_li;

 public:
  trans_node(std::vector<SMVnode *> li, NodeMtype t)
  {
    mt = t;
    ex_li = li;
  }

  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);

  std::vector<SMVnode *> get_list() { return ex_li; }
};

class invar_node : public element_node
{
  std::vector<SMVnode *> ex_li;

 public:
  invar_node(std::vector<SMVnode *> li, NodeMtype t)
  {
    mt = t;
    ex_li = li;
  }

  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
  std::vector<SMVnode *> get_list() { return ex_li; }
};

class invarspec_node : public element_node
{
  std::vector<SMVnode *> ex_li;

 public:
  invarspec_node(std::vector<SMVnode *> li, NodeMtype t)
  {
    mt = t;
    ex_li = li;
  }

  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);

  std::vector<SMVnode *> get_list() { return ex_li; }
};

class constant : public SMVnode
{
  string in;

 public:
  constant(std::string input) { in = input; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};

class identifier : public SMVnode
{
  string in;

 public:
  identifier(std::string input) { in = input; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class par_expr : public SMVnode
{
  pono::SMVnode * ex;

 public:
  par_expr(pono::SMVnode * e) { ex = e; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class not_expr : public SMVnode
{
  pono::SMVnode * ex;

 public:
  not_expr(pono::SMVnode * e) { ex = e; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class and_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  and_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class or_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  or_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class xor_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  xor_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class xnor_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  xnor_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class imp_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  imp_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class iff_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  iff_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class eq_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  eq_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class neq_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  neq_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class lt_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  lt_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class gt_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  gt_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};

class lte_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  lte_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class gte_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  gte_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class uminus_expr : public SMVnode
{
  pono::SMVnode * ex;

 public:
  uminus_expr(pono::SMVnode * e) { ex = e; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class add_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  add_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class sub_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  sub_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class mul_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  mul_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class div_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  div_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class mod_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  mod_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class sr_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  sr_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class sl_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  sl_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};

class subscript_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  subscript_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class sel_expr : public SMVnode
{
  pono::SMVnode * ex1;
  std::string start;
  std::string end;

 public:
  sel_expr(pono::SMVnode * e1, std::string s1, std::string s2)
  {
    ex1 = e1;
    start = s1;
    end = s2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class con_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  con_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};

class read_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  read_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};

class write_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;
  pono::SMVnode * ex3;

 public:
  write_expr(pono::SMVnode * e1, pono::SMVnode * e2, pono::SMVnode * e3)
  {
    ex1 = e1;
    ex2 = e2;
    ex3 = e3;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};

class apply_expr : public SMVnode
{
  std::string fun;
  std::vector<pono::SMVnode *> args;

 public:
  apply_expr(const std::string & f, std::vector<pono::SMVnode *> a)
  {
    fun = f;
    args = a;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};

class word1_expr : public SMVnode
{
  pono::SMVnode * ex;

 public:
  word1_expr(pono::SMVnode * e) { ex = e; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class bool_expr : public SMVnode
{
  pono::SMVnode * ex;

 public:
  bool_expr(pono::SMVnode * e) { ex = e; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class toint_expr : public SMVnode
{
  pono::SMVnode * ex;

 public:
  toint_expr(pono::SMVnode * e) { ex = e; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class signed_expr : public SMVnode
{
  pono::SMVnode * ex;

 public:
  signed_expr(pono::SMVnode * e) { ex = e; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class unsigned_expr : public SMVnode
{
  pono::SMVnode * ex;

 public:
  unsigned_expr(pono::SMVnode * e) { ex = e; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class extend_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  extend_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class resize_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  resize_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class union_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  union_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class set_expr : public SMVnode
{
  pono::SMVnode * ex;

 public:
  set_expr(pono::SMVnode * e) { ex = e; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class in_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  in_expr(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class ite_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;
  pono::SMVnode * ex3;

 public:
  ite_expr(pono::SMVnode * e1, pono::SMVnode * e2, pono::SMVnode * e3)
  {
    ex1 = e1;
    ex2 = e2;
    ex3 = e3;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class floor_expr : public SMVnode
{
  pono::SMVnode * ex;

 public:
  floor_expr(pono::SMVnode * e) { ex = e; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};

class constarray_type_expr : public SMVnode
{
  pono::SMVnode * ex;
  string id;

 public:
  constarray_type_expr( string input, pono::SMVnode * e){ id = input; ex = e; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};

class constarray_word_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;
  int size;

 public:
  constarray_word_expr(int sizet, pono::type_node * e1, pono::SMVnode * e2) { size = sizet; ex1 = e1; ex2 = e2; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};

class constarray_int_expr : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  constarray_int_expr(pono::type_node * e1, pono::SMVnode * e2) { ex1 = e1; ex2 = e2; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};


class case_expr : public SMVnode
{
  std::vector<pono::SMVnode *> ex_l;

 public:
  case_expr(std::vector<pono::SMVnode *> el) { ex_l = el; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};

class case_body_ex : public SMVnode
{
  pono::SMVnode * ex1;
  pono::SMVnode * ex2;

 public:
  case_body_ex(pono::SMVnode * e1, pono::SMVnode * e2)
  {
    ex1 = e1;
    ex2 = e2;
  }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};
class next_expr : public SMVnode
{
  pono::SMVnode * ex;

 public:
  next_expr(pono::SMVnode * e) { ex = e; }
  void generate_ostream(std::string name,
                        std::string prefix,
                        std::unordered_map<string, module_node *> module_list,
                        std::unordered_map<string, string> new_prefix,
                        ostream & s);
};

}  // namespace pono
