/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Hongce Zhang
 ** This file is part of the pono project.
 ** Copyright (c) 2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief
 **
 **
 **/

#include "vcd_witness_printer.h"

#include <algorithm>
#include <fstream>
#include <iostream>

#include "frontends/btor2_encoder.h"
#include "utils/logger.h"
#include "utils/str_util.h"

using namespace smt;
using namespace std;

namespace pono {

static const char date_time_format[] = "%A %Y/%m/%d  %H:%M:%S";
// The format of header :
// $date
// %date%
// $end
// $version PONO $end
// $timescale 1 ns $end

// ------------- HELPER FUNCTIONS ------------------ //

static std::vector<std::string> split(const std::string & str,
                                      const std::string & delim)
{
  std::vector<std::string> tokens;
  size_t prev = 0, pos = 0;
  do {
    pos = str.find(delim, prev);
    if (pos == std::string::npos) pos = str.length();
    std::string token = str.substr(prev, pos - prev);
    if (!token.empty()) tokens.push_back(token);
    prev = pos + delim.length();
  } while (pos < str.length() && prev < str.length());
  return tokens;
}

// convert a widith to a verilog string
static std::string width2range(uint64_t w)
{
  if (w > 1) return std::string("[") + std::to_string(w - 1) + ":0]";
  return "";
}

// copied from btor2_witness_printer,
// so that we don't need to include
// because its header contains also
// implementation, this creates troubles
static std::string as_bits(std::string val)
{
  // TODO: this makes assumptions on format of value from boolector
  //       to support other solvers, we need to be more general
  std::string res = val;

  if (val.length() < 2) {
    throw PonoException("Don't know how to interpret value: " + val);
  }

  if (res.substr(0, 2) == "#b") {
    // #b prefix -> b
    res = res.substr(1, val.length() - 1);
  } else if (res.substr(0, 2) == "#x") {
    throw PonoException("Not supporting hexadecimal format yet.");
  } else {
    res = res.substr(5, res.length() - 5);
    std::istringstream iss(res);
    std::vector<std::string> tokens(std::istream_iterator<std::string>{ iss },
                                    std::istream_iterator<std::string>());

    if (tokens.size() != 2) {
      throw PonoException("Failed to interpret " + val);
    }

    res = tokens[0];
    // get rid of ")"
    std::string width_str = tokens[1].substr(0, tokens[1].length() - 1);
    size_t width = std::stoull(width_str);
    mpz_class cval(res);
    res = cval.get_str(2);
    size_t strlen = res.length();

    if (strlen < width) {
      // pad with zeros
      res = std::string(width - strlen, '0') + res;
    } else if (strlen > width) {
      // remove prepended zeros
      res = res.erase(0, strlen - width);
    }
    return res;
  }
  return res;
}

// convert boolector value to decimal
static std::string as_decimal(std::string val)
{
  // TODO: this makes assumptions on format of value from boolector
  //       to support other solvers, we need to be more general
  std::string res = val;

  if (val.length() < 2) {
    throw PonoException("Don't know how to interpret value: " + val);
  }

  if (res.substr(0, 2) == "#b") {
    // #b prefix -> b
    res = res.substr(2, val.length() - 1);
    mpz_class cval(res, 2);
    res = cval.get_str(10);
  } else if (res.substr(0, 2) == "#x") {
    throw PonoException("Not supporting hexadecimal format yet.");
  } else {
    res = res.substr(5, res.length() - 5);
    std::istringstream iss(res);
    std::vector<std::string> tokens(std::istream_iterator<std::string>{ iss },
                                    std::istream_iterator<std::string>());

    if (tokens.size() != 2) {
      throw PonoException("Failed to interpret " + val);
    }

    res = tokens[0];
    // get rid of ")"
    std::string width_str = tokens[1].substr(0, tokens[1].length() - 1);
    size_t width = std::stoull(width_str);
    mpz_class cval(res);
    res = cval.get_str(10);
    return res;
  }
  return res;
}

// ------------- CLASS FUNCTIONS ------------------ //

VCDWitnessPrinter::VCDWitnessPrinter(
    const TransitionSystem & ts,
    const std::vector<smt::UnorderedTermMap> & cex,
    const std::unordered_map<std::string, std::string> & symbol_map)
    : inputs_(ts.inputvars()),
      states_(ts.statevars()),
      named_terms_(ts.named_terms()),
      cex_(cex),
      hash_id_cnt_(0),
      property_id_cnt_(0)
{
  // figure out the variables and their scopes
  for (auto && name_term_pair : named_terms_) {
    // It seems that next_var now will also go into named_terms
    if (ts.is_next_var(name_term_pair.second)) continue;

    bool is_reg =
        std::find(states_.begin(), states_.end(), name_term_pair.second)
        != states_.end();
    auto sk = name_term_pair.second->get_sort()->get_sort_kind();
    if (sk == smt::ARRAY) {
      continue;  // let's not worry about array so far
    }
    check_insert_scope(lookup_or_key(symbol_map, name_term_pair.first),
                       is_reg,
                       name_term_pair.second);
  }

  for (auto && state : states_) {
    if (state->get_sort()->get_sort_kind() == smt::ARRAY) {
      // we need to preprocess cex to know the indices
      std::unordered_set<std::string> indices;
      bool has_default_value = false;

      for (auto && valmap : cex) {  // for each step
        auto array_assign_pos = valmap.find(state);
        if (array_assign_pos == valmap.end())
          continue;  // we find no assignment at this step
        // should not be the case, but we skip anyway
        // peel the (store (store ...) ), find the indices

        smt::Term tmp = array_assign_pos->second;
        smt::TermVec store_children(3);
        while (tmp->get_op() == smt::Store) {
          int num = 0;
          for (auto c : tmp) {
            store_children[num] = c;
            num++;
          }
          indices.insert(as_decimal(store_children[1]->to_string()));
          tmp = store_children[0];
        }

        if (tmp->get_op().is_null() && tmp->is_value())
          has_default_value = true;
      }  // for each frame, collect
      check_insert_scope_array(lookup_or_key(symbol_map, state->to_string()),
                               indices,
                               has_default_value,
                               state);
    } else
      check_insert_scope(
          lookup_or_key(symbol_map, state->to_string()), true, state);
  }

  for (auto && input : inputs_) {
    if (input->get_sort()->get_sort_kind() == smt::ARRAY) continue;
    check_insert_scope(
        lookup_or_key(symbol_map, input->to_string()), false, input);
  }

}  // VCDWitnessPrinter -- constructor

void VCDWitnessPrinter::debug_dump() const
{
  for (uint64_t fidx = 0; fidx < cex_.size(); ++fidx) {
    logger.log(3, "------------- CEX : F{} -----------------", fidx);
    for (auto && t : cex_.at(fidx)) {
      logger.log(3, "{} -> {}", t.first->to_string(), t.second->to_string());
    }
  }
}  // debug_dump

std::string VCDWitnessPrinter::new_hash_id()
{
  return "v" + std::to_string(hash_id_cnt_++);
}

std::string VCDWitnessPrinter::new_property_id()
{
  return "assert(property" + std::to_string(property_id_cnt_++) + ")";
}

bool static is_bad_state_pattern(const std::string & n)
{
  auto dot_pos = n.rfind(':');
  for (auto pos = dot_pos + 1; pos < n.length(); ++pos) {
    char ch = n.at(pos);
    if (isdigit(ch) || ch == '.' || ch == '-') continue;
    return false;
  }
  return true;
}

void VCDWitnessPrinter::check_insert_scope(std::string full_name,
                                           bool is_reg,
                                           const smt::Term & ast)
{
  // yosys could use $... for internal unnamed nodes,
  // which maybe we don't want at all

  if (full_name.front() == '$') return;
  if (is_bad_state_pattern(full_name)) full_name = new_property_id();

  // clang-format off
  // yosys use " ; " as the separator for comment
  // HZ: I'm actually surprised that text after ' ; '
  // is not parsed by Btor2 frontend
  // so the check is mostly unuseful.
  // But it does not hurt to have it.
  // There is just one exception, that is the
  // name after bad state:
  // Yosys will output something like this:
  //
  //   155 bad 154 ./ridecore-src-buggy/topsim.v:101.13-112.8|./ridecore-src-buggy/pipeline.v:2005.11-2006.28
  //
  // and the symbols and dots will overwhelm the later code that
  // tries to sort out the hierarchy of the signal.
  // Actually this is not signal name at all.
  // That's why I use `new_property_id` above to replace it
  // clang-format on

  auto pos = full_name.find(" ; ");
  if (pos != full_name.npos) full_name = full_name.substr(0, pos);
  // gtkwave doesn't like colons in name
  std::replace(full_name.begin(), full_name.end(), ':', '_');
  auto scopes = split(full_name, ".");
  VCDScope * root = &root_scope_;
  for (size_t idx = 0; idx < scopes.size() - 1; ++idx) {
    const auto & next_scope = scopes.at(idx);
    auto pos = root->subscopes.find(next_scope);
    if (pos != root->subscopes.end()) {  // we find it
      root = &(pos->second);
    } else {  // we need to insert this scope
      root->subscopes.emplace(next_scope, VCDScope());
      root = &(root->subscopes.at(next_scope));
    }
  }  // at the end of this loop, we are at the scope to insert our variable
  const auto & short_name = scopes.back();
  uint64_t width = ast->get_sort()->get_width();

  std::map<std::string, VCDSignal> & signal_set =
      is_reg ? root->regs : root->wires;

  if (signal_set.find(short_name) != signal_set.end()) {
    // this can happen if the term is registered both under
    // named_terms_ and statevars/inputvars
    // we can ignore this issue.
    return;
  }
  auto hashid = new_hash_id();
  signal_set.emplace(
      short_name,
      VCDSignal(
          short_name + width2range(width), full_name, hashid, ast, width));
  allsig_bv_.push_back(&(signal_set.at(short_name)));
}  // end of check_insert_scope

void VCDWitnessPrinter::check_insert_scope_array(
    std::string full_name,
    const std::unordered_set<std::string> & indices,
    bool has_default,
    const smt::Term & ast)
{
  // vcd doesn't like colons in name
  std::replace(full_name.begin(), full_name.end(), ':', '_');
  auto scopes = split(full_name, ".");
  VCDScope * root = &root_scope_;
  for (size_t idx = 0; idx < scopes.size() - 1; ++idx) {
    const auto & next_scope = scopes.at(idx);
    auto pos = root->subscopes.find(next_scope);
    if (pos != root->subscopes.end()) {  // we find it
      root = &(pos->second);
    } else {  // we need to insert this scope
      root->subscopes.emplace(next_scope, VCDScope());
      root = &(root->subscopes.at(next_scope));
    }
  }  // at the end of this loop, we are at the scope to insert our variable
  const auto & short_name = scopes.back();
  uint64_t data_width = ast->get_sort()->get_elemsort()->get_width();

  std::map<std::string, VCDArray> & signal_set = root->arrays;

  if (signal_set.find(short_name) != signal_set.end()) {
    // I actually maybe should use `assert(false)` here, because
    // the code now should guarantee this will not happen
    throw PonoException(full_name + " has been registered already");
  }

  signal_set.emplace(short_name,
                     VCDArray(short_name, full_name, ast, data_width));
  auto & indices2hash = signal_set.at(short_name).indices2hash;
  for (const auto & index : indices) indices2hash.emplace(index, new_hash_id());
  if (has_default) indices2hash.emplace("default", new_hash_id());
  allsig_array_.push_back(&(signal_set.at(short_name)));
  // to do: add indices and their hashes
}  // end of check_insert_scope_array

void VCDWitnessPrinter::dump_current_scope(std::ostream & fout,
                                           const VCDScope * scope) const
{
  for (auto && r : scope->regs) {
    fout << "$var reg " << r.second.data_width << " " << r.second.hash << " "
         << r.second.vcd_name << " $end" << std::endl;
  }
  for (auto && w : scope->wires) {
    fout << "$var wire " << w.second.data_width << " " << w.second.hash << " "
         << w.second.vcd_name << " $end" << std::endl;
  }
  for (auto && a : scope->arrays) {
    auto data_width = a.second.data_width;
    for (auto && idx_hash_pair : a.second.indices2hash)
      fout << "$var reg " << data_width << " " << idx_hash_pair.second << " "
           << a.second.vcd_name + "[" + idx_hash_pair.first + "]"
                  + width2range(data_width)
           << " $end" << std::endl;
  }
  // let's go for the submodules
  for (auto pos = scope->subscopes.begin(); pos != scope->subscopes.end();
       ++pos) {
    fout << "$scope module " << pos->first << " $end" << std::endl;
    dump_current_scope(fout, &(pos->second));
    fout << "$upscope $end" << std::endl;
  }
}  // end of dump_current_scope

void VCDWitnessPrinter::DumpScopes(std::ostream & fout) const
{
  dump_current_scope(fout, &root_scope_);
}

void VCDWitnessPrinter::GenHeader(std::ostream & fout) const
{
  fout << "$date" << std::endl;
  {
    char buffer[100];
    time_t rawtime;
    struct tm * timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    if (strftime(buffer, 100, date_time_format, timeinfo) == 0)
      throw PonoException("Bug: time2string conversion failed.");
    fout << buffer << std::endl;
  }
  fout << "$end" << std::endl;
  fout << "$version PONO $end" << std::endl;
  fout << "$timescale 1 ns $end" << std::endl;
  DumpScopes(fout);
  fout << "$enddefinitions $end" << std::endl;
}  // end of GenHeader

void VCDWitnessPrinter::dump_all(
    const smt::UnorderedTermMap & valmap,
    std::unordered_map<std::string, std::string> & valbuf,
    uint64_t t,
    std::ostream & fout) const
{
  for (auto && sig_bv_ptr : allsig_bv_) {
    auto pos = valmap.find(sig_bv_ptr->ast);
    if (pos == valmap.end()) {
      logger.log(1,
                 "missing value in provided trace @{}: {}",
                 t,
                 sig_bv_ptr->full_name);
      continue;
    }
    auto val = as_bits(pos->second->to_string());
    valbuf.emplace(sig_bv_ptr->hash, val);
    fout << val << " " << sig_bv_ptr->hash << std::endl;
  }  // for all bv signals

  for (auto && sig_array_ptr : allsig_array_) {
    auto pos = valmap.find(sig_array_ptr->ast);
    if (pos == valmap.end()) {
      logger.log(1,
                 "missing value in provided trace @{}: {}",
                 t,
                 sig_array_ptr->full_name);
      continue;
    }
    smt::Term memvalue = pos->second;
    smt::TermVec store_children(3);
    while (memvalue->get_op() == smt::Store) {
      int num = 0;
      for (auto c : memvalue) {
        store_children[num] = c;
        num++;
      }

      auto addr = as_decimal(store_children[1]->to_string());
      auto data = as_bits(store_children[2]->to_string());
      auto addr_pos = sig_array_ptr->indices2hash.find(addr);
      if (addr_pos != sig_array_ptr->indices2hash.end()) {
        valbuf.emplace(addr_pos->second, data);
        fout << data << " " << addr_pos->second << std::endl;
      } else {
        logger.log(1,
                   "missing addr index for array: {}: , addr : {}",
                   sig_array_ptr->full_name,
                   addr);
      }

      memvalue = store_children[0];
    }

    if (memvalue->get_op().is_null() && memvalue->is_value()) {
      smt::Term const_val = *(memvalue->begin());
      auto data_default = as_bits(const_val->to_string());
      auto addr_pos = sig_array_ptr->indices2hash.find("default");
      if (addr_pos != sig_array_ptr->indices2hash.end()) {
        valbuf.emplace(addr_pos->second, data_default);
        fout << data_default << " " << addr_pos->second << std::endl;
      } else {
        logger.log(1,
                   "missing addr index for array: {}: , addr : {}",
                   sig_array_ptr->full_name,
                   "-default-");
      }
    }  // handling the inner constant default
  }  // for all array signals
  // TODO : mems
}  // end of VCDWitnessPrinter::dump_all

void VCDWitnessPrinter::dump_diff(
    const smt::UnorderedTermMap & valmap,
    std::unordered_map<std::string, std::string> & valprev,
    uint64_t t,
    std::ostream & fout) const
{
  for (auto && sig_bv_ptr : allsig_bv_) {
    auto pos = valmap.find(sig_bv_ptr->ast);
    if (pos == valmap.end()) {
      logger.log(1,
                 "missing value in provided trace @{}: {}",
                 t,
                 sig_bv_ptr->full_name);
      continue;
    }
    auto val = as_bits(pos->second->to_string());
    auto prev_pos = valprev.find(sig_bv_ptr->hash);
    if (prev_pos == valprev.end()) {
      valprev.emplace(sig_bv_ptr->hash, val);
      fout << val << " " << sig_bv_ptr->hash << std::endl;
      logger.log(1,
                 "Bug, {} was not cached before time : {}.",
                 sig_bv_ptr->full_name,
                 std::to_string(t));
      continue;
    }
    if (prev_pos->second == val) continue;
    // update old value and print
    prev_pos->second = val;
    fout << val << " " << sig_bv_ptr->hash << std::endl;
  }  // for all bv signals

  for (auto && sig_array_ptr : allsig_array_) {
    auto pos = valmap.find(sig_array_ptr->ast);
    if (pos == valmap.end()) {
      logger.log(1,
                 "missing value in provided trace @{}: {}",
                 t,
                 sig_array_ptr->full_name);
      continue;
    }
    smt::Term memvalue = pos->second;
    smt::TermVec store_children(3);
    while (memvalue->get_op() == smt::Store) {  // peel the (store (store ...))
      int num = 0;
      for (auto c : memvalue) {
        store_children[num] = c;
        num++;
      }

      auto addr = as_decimal(store_children[1]->to_string());
      auto data = as_bits(store_children[2]->to_string());
      auto addr_pos = sig_array_ptr->indices2hash.find(addr);
      if (addr_pos != sig_array_ptr->indices2hash.end()) {
        auto prev_pos = valprev.find(addr_pos->second);
        if (prev_pos == valprev.end()) {
          valprev.emplace(addr_pos->second, data);
          fout << data << " " << addr_pos->second << std::endl;
          // this happens if some elements are not assigned by
          // the solver in the beginning, I think this is also
          // common and GtkWave is able to handle it (by having X, don't care)
          // before it is first assigned
          // so you can consider even remove the logger below
          logger.log(3,
                     "{} was not cached before time : {}.",
                     sig_array_ptr->full_name + "[" + addr + "]",
                     std::to_string(t));
        } else {
          if (prev_pos->second != data) {
            prev_pos->second = data;  // update the value
            fout << data << " " << addr_pos->second << std::endl;
          }
        }  // exists in prev pos or not
      } else {
        // we should have already recorded all the indices by
        // going through all the frames earlier on
        // so if the case happens, it should be a bug actually
        logger.log(1,
                   "missing addr index for array: {}: , addr : {}",
                   sig_array_ptr->full_name,
                   addr);
      }
      memvalue = store_children[0];
    }

    if (memvalue->get_op().is_null() && memvalue->is_value()) {
      smt::Term const_val = *(memvalue->begin());
      auto data_default = as_bits(const_val->to_string());
      auto addr_pos = sig_array_ptr->indices2hash.find("default");
      if (addr_pos != sig_array_ptr->indices2hash.end()) {
        auto prev_pos = valprev.find(addr_pos->second);
        if (prev_pos == valprev.end()) {
          valprev.emplace(addr_pos->second, data_default);
          fout << data_default << " " << addr_pos->second << std::endl;
          // this happens if some elements are not assigned by
          // the solver in the beginning, I think this is also
          // common and GtkWave is able to handle it (by having X, don't care)
          // before it is first assigned
          // so you can consider even remove the logger below
          logger.log(3,
                     "{} was not cached before time : {}.",
                     sig_array_ptr->full_name + "[default]",
                     std::to_string(t));
        } else {
          if (prev_pos->second != data_default) {
            prev_pos->second = data_default;  // update the value
            fout << data_default << " " << addr_pos->second << std::endl;
          }
        }  // exists in prev pos or not
      } else {
        // we should have already recorded all the indices by
        // going through all the frames earlier on
        // so if the case happens, it should be a bug actually
        logger.log(1,
                   "missing addr index for array: {}: , addr : {}",
                   sig_array_ptr->full_name,
                   "-default-");
      }
    }  // handling the inner constant default
  }  // for all array signals
}  // end of VCDWitnessPrinter::dump_diff

void VCDWitnessPrinter::DumpValues(std::ostream & fout) const
{
  // at time 0 we dump all the values
  // and then at each later time, we compare and
  // see if there is any difference, if not
  // we just skip it
  // finally add an empty time-tick
  if (cex_.empty()) throw PonoException("No trace to dump");

  std::unordered_map<std::string, std::string> hash_to_value_map;
  // used to store the previous value for comparison
  fout << "#0" << std::endl;
  dump_all(cex_.at(0), hash_to_value_map, 0, fout);
  for (uint64_t t = 1; t < cex_.size(); ++t) {
    fout << "#" << t << std::endl;
    dump_diff(cex_.at(t), hash_to_value_map, t, fout);
  }
  fout << "#" << cex_.size() << std::endl;
}

void VCDWitnessPrinter::dump_trace_to_file(
    const std::string & vcd_file_name) const
{
  std::ofstream fout(vcd_file_name);
  if (!fout.is_open())
    throw PonoException("Unable to write to : " + vcd_file_name);

  GenHeader(fout);
  DumpValues(fout);
  logger.log(0, "Trace written to " + vcd_file_name);
}  // dump_trace_to_file

}  // namespace pono
