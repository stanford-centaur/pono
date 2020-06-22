#include "smv_node.h"

void cosa::module_node::pre(
    std::unordered_map<std::string, module_node *> module_list,
    std::ostream & s)
{
  std::unordered_map<string, string> new_prefix;
  new_prefix[module_name] = "";
  std::string prefix = "";
  var_li->preprocess(module_name, "", module_list, &new_prefix, s);
  var_li->generate_ostream(module_name, prefix, module_list, new_prefix, s);
  ivar_li->generate_ostream(module_name, prefix, module_list, new_prefix, s);
  frozenvar_li->generate_ostream(
      module_name, prefix, module_list, new_prefix, s);
  define_li->generate_ostream(module_name, prefix, module_list, new_prefix, s);
  assign_li->generate_ostream(module_name, prefix, module_list, new_prefix, s);
  init_li->generate_ostream(module_name, prefix, module_list, new_prefix, s);
  trans_li->generate_ostream(module_name, prefix, module_list, new_prefix, s);
  invar_li->generate_ostream(module_name, prefix, module_list, new_prefix, s);
  invarspec_li->generate_ostream(
      module_name, prefix, module_list, new_prefix, s);
}

void cosa::module_node::preprocess(
    std::string parent,
    std::string prefix,
    std::unordered_map<string, string> * new_prefix,
    std::unordered_map<std::string, module_node *> module_list,
    std::vector<SMVnode *> id_li,
    std::ostream & s)
{
  if (par_li.size() != id_li.size()) {
    throw CosaException("module declaration wrong");
  }
  for (int i = 0; i < id_li.size(); i++) {
    std::string s = par_li[i];
    new_par[s] = id_li[i];
  }
  (*new_prefix)[module_name] = prefix;
  par_name = parent;
  std::unordered_map<string, string> prefix_li = *new_prefix;
  var_li->preprocess(module_name, prefix, module_list, new_prefix, s);
  var_li->generate_ostream(module_name, prefix, module_list, prefix_li, s);
  ivar_li->generate_ostream(module_name, prefix, module_list, prefix_li, s);
  frozenvar_li->generate_ostream(
      module_name, prefix, module_list, prefix_li, s);
  define_li->generate_ostream(module_name, prefix, module_list, prefix_li, s);
  assign_li->generate_ostream(module_name, prefix, module_list, prefix_li, s);
  init_li->generate_ostream(module_name, prefix, module_list, prefix_li, s);
  trans_li->generate_ostream(module_name, prefix, module_list, prefix_li, s);
  invar_li->generate_ostream(module_name, prefix, module_list, prefix_li, s);
  invarspec_li->generate_ostream(
      module_name, prefix, module_list, prefix_li, s);
}

void cosa::type_node::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << type_name;
}

void cosa::var_node_c::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  if (vt != ModuleT) {
    unordered_map<string, SMVnode *> new_par =
        module_list[name]->get_namelist();
    if (new_par.find(id) != new_par.end())
      throw CosaException("duplicately defined");
    s << new_prefix[name] << id << ":" << type << ";" << endl;
  }
}

void cosa::ivar_node_c::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  unordered_map<string, SMVnode *> new_par = module_list[name]->get_namelist();
  if (new_par.find(id) != new_par.end())
    throw CosaException("duplicately defined");
  s << new_prefix[name] << id << " : " << type << " ; " << endl;
}

void cosa::frozenvar_node_c::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  unordered_map<string, SMVnode *> new_par = module_list[name]->get_namelist();
  if (new_par.find(id) != new_par.end())
    throw CosaException("duplicately defined");
  s << new_prefix[name] << id << " : " << type << " ; " << endl;
}

void cosa::define_node_c::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  unordered_map<string, SMVnode *> new_par = module_list[name]->get_namelist();
  if (new_par.find(id) != new_par.end())
    throw CosaException("duplicately defined");
  s << prefix << id << " := ";
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ;" << endl;
}

void cosa::assign_node_c::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  unordered_map<string, SMVnode *> new_par = module_list[name]->get_namelist();
  if (new_par.find(id) != new_par.end())
    throw CosaException("duplicately defined");
  if (pre == "")
    s << prefix << id << " := ";
  else
    s << "( " << pre << prefix << id << ") :=";
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ;" << endl;
}

void cosa::init_node_c::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ;" << endl;
}

void cosa::trans_node_c::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ;" << endl;
}

void cosa::invar_node_c::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ;" << endl;
}

void cosa::invarspec_node_c::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ;" << endl;
}

void cosa::var_node::preprocess(
    std::string module,
    std::string prefix,
    std::unordered_map<std::string, module_node *> module_list,
    std::unordered_map<string, string> * new_prefix,
    ostream & s)
{
  for (int i = 0; i < ex_li.size(); i++) {
    var_node_c * temp = ex_li[i];
    if (temp->getVarType() == ModuleT) {
      unordered_map<string, module_node *>::iterator it =
          module_list.find(temp->getmodtype()->getName());
      if (it != module_list.end()) {
        std::string prefix = (*new_prefix)[module] + ex_li[i]->getName() + ".";
        it->second->preprocess(module,
                               prefix,
                               new_prefix,
                               module_list,
                               temp->getmodtype()->get_list(),
                               s);
      }
    }
  }
}

void cosa::var_node::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  if (!ex_li.empty()) {
    s << "VAR" << endl;
    for (int i = ex_li.size() - 1; i > -1; i--) {
      ex_li[i]->generate_ostream(name, prefix, module_list, new_prefix, s);
    }
  }
}

void cosa::ivar_node::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  if (!ex_li.empty()) {
    s << "IVAR" << endl;
    for (int i = ex_li.size() - 1; i > -1; i--) {
      ex_li[i]->generate_ostream(name, prefix, module_list, new_prefix, s);
    }
  }
}

void cosa::frozenvar_node::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  if (!ex_li.empty()) {
    s << "FROZENVAR" << endl;
    for (int i = ex_li.size() - 1; i > -1; i--) {
      ex_li[i]->generate_ostream(name, prefix, module_list, new_prefix, s);
    }
  }
}

void cosa::define_node::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  if (!ex_li.empty()) {
    s << "DEFINE" << endl;
    for (int i = ex_li.size() - 1; i > -1; i--) {
      ex_li[i]->generate_ostream(name, prefix, module_list, new_prefix, s);
    }
  }
}

void cosa::assign_node::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  if (!ex_li.empty()) {
    s << "ASSIGN" << endl;
    for (int i = ex_li.size() - 1; i > -1; i--) {
      ex_li[i]->generate_ostream(name, prefix, module_list, new_prefix, s);
    }
  }
}

void cosa::init_node::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  if (!ex_li.empty()) {
    s << "INIT" << endl;
    for (int i = ex_li.size() - 1; i > -1; i--) {
      ex_li[i]->generate_ostream(name, prefix, module_list, new_prefix, s);
    }
  }
}

void cosa::trans_node::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  if (!ex_li.empty()) {
    s << "TRANS" << endl;
    for (int i = ex_li.size() - 1; i > -1; i--) {
      ex_li[i]->generate_ostream(name, prefix, module_list, new_prefix, s);
    }
  }
}

void cosa::invar_node::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  if (!ex_li.empty()) {
    s << "INVAR" << endl;
    for (int i = ex_li.size() - 1; i > -1; i--) {
      ex_li[i]->generate_ostream(name, prefix, module_list, new_prefix, s);
    }
  }
}

void cosa::invarspec_node::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  if (!ex_li.empty()) {
    s << "INVARSPEC" << endl;
    for (int i = ex_li.size() - 1; i > -1; i--) {
      ex_li[i]->generate_ostream(name, prefix, module_list, new_prefix, s);
    }
  }
}

void cosa::constant::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << " ";
  s << in;
  s << " ";
}

void cosa::identifier::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  unordered_map<string, SMVnode *> name_list =
      module_list[name]->get_namelist();
  if (name_list.find(in) != name_list.end()) {
    std::unordered_map<string, SMVnode *>::iterator it = name_list.find(in);
    it->second->generate_ostream(
        module_list[name]->get_par(), prefix, module_list, new_prefix, s);
  } else {
    s << " " << new_prefix[name] << in << " ";
  }
}

void cosa::par_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << " ( ";
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ) ";
}

void cosa::not_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << " !";
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::and_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " & ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::or_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " | ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::xor_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " xor ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::xnor_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " xnor ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::imp_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " -> ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::iff_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " <-> ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::eq_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " = ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::neq_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " != ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::lt_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " < ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::gt_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " > ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::lte_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " <= ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::gte_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " >= ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::uminus_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << " - ";
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::add_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " + ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::sub_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " - ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::mul_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " * ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::div_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " / ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::mod_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " mod ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::sr_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " >>  ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::sl_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " << ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::subscript_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << " [ ";
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ] ";
}

void cosa::sel_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " [ ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " : ";
  ex3->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ] ";
}

void cosa::con_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " :: ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::read_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << "READ ( ";
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " , ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ) ";
}

void cosa::write_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << "WRITE ( ";
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " , ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " , ";
  ex3->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::word1_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << "word1 ( ";
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ) ";
}

void cosa::bool_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << "bool ( ";
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ) ";
}

void cosa::toint_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << "toint ( ";
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ) ";
}

void cosa::signed_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << " signed ( ";
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ) ";
}

void cosa::unsigned_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << "unsigned ( ";
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ) ";
}

void cosa::extend_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << " extend ";
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " , ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::resize_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << " resize ";
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " , ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::union_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " union ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::set_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << "!";
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::in_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " in ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::ite_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ? ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " : ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void cosa::floor_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << "floor ( ";
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ) ";
}

void cosa::case_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << " case ";
  for (int i = 0; i < ex_l.size(); i++) {
    ex_l[i]->generate_ostream(name, prefix, module_list, new_prefix, s);
    s << endl;
  }
  s << " esac ";
}

void cosa::case_body_ex::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << " ";
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " : ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ; ";
}

void cosa::next_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << "next ( ";
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ) ";
}