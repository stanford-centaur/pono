#include "smv_node.h"

void pono::module_node::process_main(
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
  fun_li->generate_ostream(module_name, prefix, module_list, new_prefix, s);
  define_li->generate_ostream(module_name, prefix, module_list, new_prefix, s);
  assign_li->generate_ostream(module_name, prefix, module_list, new_prefix, s);
  init_li->generate_ostream(module_name, prefix, module_list, new_prefix, s);
  trans_li->generate_ostream(module_name, prefix, module_list, new_prefix, s);
  invar_li->generate_ostream(module_name, prefix, module_list, new_prefix, s);
  invarspec_li->generate_ostream(
      module_name, prefix, module_list, new_prefix, s);
}

void pono::module_node::preprocess(
    std::string parent,
    std::string prefix,
    std::unordered_map<string, string> * new_prefix,
    std::unordered_map<std::string, module_node *> module_list,
    std::vector<SMVnode *> id_li,
    std::ostream & s)
{
  if (par_li.size() != id_li.size()) {
    throw PonoException("module declaration wrong");
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
  fun_li->generate_ostream(module_name, prefix, module_list, prefix_li, s);
  define_li->generate_ostream(module_name, prefix, module_list, prefix_li, s);
  assign_li->generate_ostream(module_name, prefix, module_list, prefix_li, s);
  init_li->generate_ostream(module_name, prefix, module_list, prefix_li, s);
  trans_li->generate_ostream(module_name, prefix, module_list, prefix_li, s);
  invar_li->generate_ostream(module_name, prefix, module_list, prefix_li, s);
  invarspec_li->generate_ostream(
      module_name, prefix, module_list, prefix_li, s);
}

void pono::type_node::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << type_name;
}

void pono::var_node_c::generate_ostream(
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
      throw PonoException("duplicately defined");
    s << new_prefix[name] << id << ":" << type << ";" << endl;
  }
}

void pono::ivar_node_c::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  unordered_map<string, SMVnode *> new_par = module_list[name]->get_namelist();
  if (new_par.find(id) != new_par.end())
    throw PonoException("duplicately defined");
  s << new_prefix[name] << id << " : " << type << " ; " << endl;
}

void pono::frozenvar_node_c::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  unordered_map<string, SMVnode *> new_par = module_list[name]->get_namelist();
  if (new_par.find(id) != new_par.end())
    throw PonoException("duplicately defined");
  s << new_prefix[name] << id << " : " << type << " ; " << endl;
}

void pono::fun_node_c::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  unordered_map<string, SMVnode *> new_par = module_list[name]->get_namelist();
  if (new_par.find(id) != new_par.end())
    throw PonoException("duplicately defined");
  s << new_prefix[name] << id << " : " << type << " ; " << endl;
}

void pono::define_node_c::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  unordered_map<string, SMVnode *> new_par = module_list[name]->get_namelist();
  if (new_par.find(id) != new_par.end())
    throw PonoException("duplicately defined");
  s << prefix << id << " := ";
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ;" << endl;
}

void pono::assign_node_c::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  unordered_map<string, SMVnode *> new_par = module_list[name]->get_namelist();
  if (new_par.find(id) != new_par.end())
    throw PonoException("duplicately defined");
  if (pre == "")
    s << prefix << id << " := ";
  else
    s << pre << " ( " << prefix << id << ") :=";
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ;" << endl;
}

void pono::init_node_c::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ;" << endl;
}

void pono::trans_node_c::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ;" << endl;
}

void pono::invar_node_c::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ;" << endl;
}

void pono::invarspec_node_c::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ;" << endl;
}

void pono::var_node::preprocess(
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

void pono::var_node::generate_ostream(
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

void pono::ivar_node::generate_ostream(
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

void pono::frozenvar_node::generate_ostream(
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

void pono::fun_node::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  if (!ex_li.empty()) {
    s << "FUN" << endl;
    for (int i = ex_li.size() - 1; i > -1; i--) {
      ex_li[i]->generate_ostream(name, prefix, module_list, new_prefix, s);
    }
  }
}

void pono::define_node::generate_ostream(
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

void pono::assign_node::generate_ostream(
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

void pono::init_node::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  if (!ex_li.empty()) {
    for (int i = ex_li.size() - 1; i > -1; i--) {
      s << "INIT" << endl;
      ex_li[i]->generate_ostream(name, prefix, module_list, new_prefix, s);
    }
  }
}

void pono::trans_node::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  if (!ex_li.empty()) {
    for (int i = ex_li.size() - 1; i > -1; i--) {
      s << "TRANS" << endl;
      ex_li[i]->generate_ostream(name, prefix, module_list, new_prefix, s);
    }
  }
}

void pono::invar_node::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  if (!ex_li.empty()) {
    for (int i = ex_li.size() - 1; i > -1; i--) {
      s << "INVAR" << endl;
      ex_li[i]->generate_ostream(name, prefix, module_list, new_prefix, s);
    }
  }
}

void pono::invarspec_node::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  if (!ex_li.empty()) {
    for (int i = ex_li.size() - 1; i > -1; i--) {
      s << "INVARSPEC" << endl;
      ex_li[i]->generate_ostream(name, prefix, module_list, new_prefix, s);
    }
  }
}

void pono::constant::generate_ostream(
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

void pono::identifier::generate_ostream(
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

void pono::par_expr::generate_ostream(
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

void pono::not_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << " !";
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void pono::and_expr::generate_ostream(
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

void pono::or_expr::generate_ostream(
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

void pono::xor_expr::generate_ostream(
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

void pono::xnor_expr::generate_ostream(
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

void pono::imp_expr::generate_ostream(
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

void pono::iff_expr::generate_ostream(
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

void pono::eq_expr::generate_ostream(
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

void pono::neq_expr::generate_ostream(
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

void pono::lt_expr::generate_ostream(
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

void pono::gt_expr::generate_ostream(
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

void pono::lte_expr::generate_ostream(
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

void pono::gte_expr::generate_ostream(
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

void pono::uminus_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << " - ";
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void pono::add_expr::generate_ostream(
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

void pono::sub_expr::generate_ostream(
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

void pono::mul_expr::generate_ostream(
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

void pono::div_expr::generate_ostream(
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

void pono::mod_expr::generate_ostream(
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

void pono::sr_expr::generate_ostream(
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

void pono::sl_expr::generate_ostream(
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

void pono::subscript_expr::generate_ostream(
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

void pono::sel_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " [ " << start;
  s << " : " << end;
  s << " ] ";
}

void pono::con_expr::generate_ostream(
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

void pono::read_expr::generate_ostream(
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

void pono::write_expr::generate_ostream(
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
  s << " ) ";
}

void pono::apply_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << fun << "(";
  size_t argsize = args.size();
  for (size_t i = 0; i + 1 < argsize; ++i) {
    args[i]->generate_ostream(name, prefix, module_list, new_prefix, s);
    s << ", ";
  }
  args[argsize - 1]->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << ")";
}

void pono::word1_expr::generate_ostream(
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

void pono::bool_expr::generate_ostream(
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

void pono::toint_expr::generate_ostream(
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

void pono::signed_expr::generate_ostream(
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

void pono::unsigned_expr::generate_ostream(
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

void pono::extend_expr::generate_ostream(
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

void pono::resize_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << " resize (";
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " , ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ) ";
}

void pono::union_expr::generate_ostream(
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

void pono::set_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << "!";
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void pono::in_expr::generate_ostream(
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

void pono::ite_expr::generate_ostream(
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
  ex3->generate_ostream(name, prefix, module_list, new_prefix, s);
}

void pono::floor_expr::generate_ostream(
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

void pono::constarray_type_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << "CONSTARRAY ( typeof ( " << prefix << id << " ) ,";
  ex->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ) ";
}

void pono::constarray_word_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << "CONSTARRAY ( array word " << size << " of";
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " , ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ) ";
}

void pono::constarray_int_expr::generate_ostream(
    std::string name,
    std::string prefix,
    std::unordered_map<string, module_node *> module_list,
    std::unordered_map<string, string> new_prefix,
    ostream & s)
{
  s << "CONSTARRAY ( array integer of ";
  ex1->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " , ";
  ex2->generate_ostream(name, prefix, module_list, new_prefix, s);
  s << " ) ";
}

void pono::case_expr::generate_ostream(
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

void pono::case_body_ex::generate_ostream(
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

void pono::next_expr::generate_ostream(
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
