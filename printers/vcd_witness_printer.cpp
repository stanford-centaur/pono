/*********************                                                        */
/*! \file 
 ** \verbatim
 ** Top contributors (to current version):
 **   Hongce Zhang
 ** This file is part of the cosa2 project.
 ** Copyright (c) 2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief 
 **
 ** 
 **/

#include "utils/logger.h"
#include "frontends/btor2_encoder.h"
#include "smt-switch/boolector_factory.h"

#include "vcd_witness_printer.h"
// #include "btor2_witness_printer.h"

#include <iostream>
#include <fstream>

using namespace smt;
using namespace std;

namespace cosa {

static const char date_time_format [] = "%A %Y/%m/%d  %H:%M:%S";
// const std::string_view vcd_header (R"**(
// $date
// %date%
// $end
// $version COSA2 $end
// $timescale 1 ns $end
// )**");

// ------------- HELPER FUNCTIONS ------------------ //

static std::vector<std::string> split(const std::string& str,
                               const std::string& delim) {
  std::vector<std::string> tokens;
  size_t prev = 0, pos = 0;
  do {
    pos = str.find(delim, prev);
    if (pos == std::string::npos)
      pos = str.length();
    std::string token = str.substr(prev, pos - prev);
    if (!token.empty())
      tokens.push_back(token);
    prev = pos + delim.length();
  } while (pos < str.length() && prev < str.length());
  return tokens;
}


/// convert a widith to a verilog string
static std::string width2range(uint64_t w) {
  if (w > 1)
    return std::string("[") + std::to_string(w - 1) + ":0]";
  return "";
}

/// copied from btor2_witness_printer,
/// so that we don't need to include
/// because its header contains also 
/// implementation, this creates troubles
static std::string as_bits(std::string val)
{
  // TODO: this makes assumptions on format of value from boolector
  //       to support other solvers, we need to be more general
  std::string res = val;

  if (val.length() < 2) {
    throw CosaException("Don't know how to interpret value: " + val);
  }

  if (res.substr(0, 2) == "#b") {
    // #b prefix -> b
    res = res.substr(1, val.length() - 1);
  } else if (res.substr(0, 2) == "#x") {
    throw CosaException("Not supporting hexadecimal format yet.");
  } else {
    res = res.substr(5, res.length() - 5);
    std::istringstream iss(res);
    std::vector<std::string> tokens(std::istream_iterator<std::string>{ iss },
                                    std::istream_iterator<std::string>());

    if (tokens.size() != 2) {
      throw CosaException("Failed to interpret " + val);
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


// ------------- CLASS FUNCTIONS ------------------ //


VCDWitnessPrinter::VCDWitnessPrinter(const BTOR2Encoder & btor_enc,
                    const TransitionSystem & ts) :
  inputs_(btor_enc.inputsvec()),
  states_(btor_enc.statesvec()),
  no_next_states_(btor_enc.no_next_states()),
  has_states_without_next_(!no_next_states_.empty()),
  named_terms_(ts.named_terms()),
  hash_id_cnt_(0)
{
  // figure out the variables and their scopes
  for (auto && name_term_pair : named_terms_) {
    bool is_reg = 
      std::find(states_.begin(), states_.end(), name_term_pair.second)
      != states_.end();
    auto sk = name_term_pair.second->get_sort()->get_sort_kind();
    if (sk == smt::ARRAY) {
      continue; // let's not worry about array so far
    }
    check_insert_scope(name_term_pair.first, is_reg, name_term_pair.second);
  }

  for (auto && state : states_) {
    if(state->get_sort()->get_sort_kind() == smt::ARRAY)
      continue;
    check_insert_scope(state->to_string(), true, state);
  }


  for (auto && input : inputs_) {
    if(input->get_sort()->get_sort_kind() == smt::ARRAY)
      continue;
    check_insert_scope(input->to_string(), false, input);
  }


  // you also need to figure out the indices used for each array
  // and dump these values when needed and dump the default value
  logger.log(3, "-------------Input Dump-----------------");
  for (auto && i: inputs_) {
    logger.log(3, "{}", i->to_string());
  }
  logger.log(3, "-------------State Dump-----------------");
  for (auto && s: states_) {
    logger.log(3, "{}", s->to_string());
  }
  logger.log(3, "-------------No Next States-----------------");
  for (auto && s: no_next_states_) {
    logger.log(3, "{}", s.second->to_string());
  }
}

void VCDWitnessPrinter::DebugDump(const std::vector<smt::UnorderedTermMap> & cex) const {
  for (uint64_t fidx = 0; fidx < cex.size(); ++ fidx) {
    logger.log(3, "------------- CEX : F{} -----------------", fidx);
    for (auto && t : cex.at(fidx)) {
      logger.log(3, "{} -> {}", t.first->to_string(), t.second->to_string() );
    }
  }
}

std::string VCDWitnessPrinter::new_hash_id() {
  return "v" + std::to_string(hash_id_cnt_++);
}

void VCDWitnessPrinter::check_insert_scope(const std::string& full_name, bool is_reg,
  const smt::Term & ast)
{
  auto scopes = split(full_name, ".");
  VCDScope * root = & root_scope_;
  for (size_t idx = 0; idx < scopes.size() - 1; ++idx) {
    const auto & next_scope = scopes.at(idx);
    auto pos = root->subscopes.find(next_scope);
    if (pos != root->subscopes.end()) { // we find it
      root = & (pos->second);
    } else { // we need to insert this scope
      root->subscopes.insert(std::make_pair(next_scope, VCDScope()));
      root = &( root->subscopes.at(next_scope) );
    }
  } // at the end of this loop, we are at the scope to insert our variable
  const auto & short_name = scopes.back();
  uint64_t width = ast->get_sort()->get_width();

  std::map<std::string, VCDSignal> & signal_set = is_reg ? root->regs : root->wires;

  if (signal_set.find(short_name) != signal_set.end()) {
    // TODO: only check in Debug
    throw CosaException(full_name + " has been registered already");
  }
  auto hashid = new_hash_id();
  signal_set.insert(std::make_pair(short_name,
    VCDSignal(
      short_name + width2range(width), 
      full_name,  hashid , ast, width)));
  hash2sig_bv_.insert(std::make_pair(hashid, &(signal_set.at(short_name)) ));
}


void VCDWitnessPrinter::dump_current_scope(std::ostream & fout, const VCDScope *scope) const {
  for (auto && r : scope->regs) {
    fout << "$var reg " << r.second.data_width << " " << r.second.hash << " "
         << r.second.vcd_name << " $end" << std::endl;
  }
  for (auto && w : scope->wires) {
    fout << "$var wire " << w.second.data_width << " " << w.second.hash << " "
         << w.second.vcd_name << " $end" << std::endl;
  }
  // let's go for the submodules
  for (auto pos = scope->subscopes.begin() ; pos != scope->subscopes.end() ; ++ pos) {
    fout << "$scope module " << pos->first << " $end" << std::endl;
    dump_current_scope(fout, &(pos->second));
    fout << "$upscope $end" << std::endl;
  }
} // end of dump_current_scope

void VCDWitnessPrinter::DumpScopes(std::ostream & fout) const {
  dump_current_scope(fout, &root_scope_);
}

void VCDWitnessPrinter::GenHeader(std::ostream & fout) const {
  fout << "$date" << std::endl;
  {
    char buffer [100];
    time_t rawtime;
    struct tm * timeinfo;
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    if ( strftime(buffer, 100, date_time_format, timeinfo) == 0)
      throw CosaException("Bug: time2string conversion failed.");
    fout << buffer << std::endl;
  }
  fout << "$end" << std::endl;
  fout << "$version CoSA2 $end" << std::endl;
  fout << "$timescale 1 ns $end" << std::endl;
  DumpScopes(fout);
  fout << "$enddefinitions $end" << std::endl;
} // GenHeader

void VCDWitnessPrinter::dump_all(const smt::UnorderedTermMap & valmap,
  std::unordered_map<std::string, std::string> & valbuf,
  uint64_t t, std::ostream & fout) const {
  for (auto && hash_sig_pair : hash2sig_bv_) {
    auto pos = valmap.find(hash_sig_pair.second->ast);
    if (pos == valmap.end()) {
      logger.log(0, "missing value in provided trace @{}: {} ,{}, {}" ,
        t,
        hash_sig_pair.second->full_name,
        hash_sig_pair.first, 
        hash_sig_pair.second->ast->to_string());
      continue;
    }
    auto val = as_bits(pos->second->to_string());
    valbuf.insert(std::make_pair(
      hash_sig_pair.first,
      val
    ));
    fout << val << " " << hash_sig_pair.first << std::endl;
  }
  // TODO : mems
}

void VCDWitnessPrinter::dump_diff(const smt::UnorderedTermMap & valmap,
  std::unordered_map<std::string, std::string> & valprev,
  uint64_t t, std::ostream & fout) const {

  for (auto && hash_sig_pair : hash2sig_bv_) {
    auto pos = valmap.find(hash_sig_pair.second->ast);
    if (pos == valmap.end()) {
      logger.log(0, "missing value in provided trace @{}: {} ,{}, {}" ,
        t,
        hash_sig_pair.second->full_name,
        hash_sig_pair.first, 
        hash_sig_pair.second->ast->to_string());
      continue;
    }
    auto val = as_bits(pos->second->to_string());
    auto prev_pos = valprev.find(hash_sig_pair.first);
    if (prev_pos == valprev.end())
      throw CosaException("Bug, "+ hash_sig_pair.first + " is not cached.");
    if ( prev_pos->second == val )
      continue;
    // update old value and print
    prev_pos->second = val;
    fout << val << " " << hash_sig_pair.first << std::endl;
  }
}

void VCDWitnessPrinter::DumpValues(std::ostream & fout, 
  const std::vector<smt::UnorderedTermMap> & cex) const {
  // at time 0 we dump all the values
  // and then at each later time, we compare and
  // see if there is any difference, if not
  // we just skip it
  // finally add an empty time-tick
  if (cex.empty())
    throw CosaException("No trace to dump");

  std::unordered_map<std::string, std::string> hash_to_value_map;
  // used to store the previous value for comparison
  fout << "#0" << std::endl;
  dump_all(cex.at(0), hash_to_value_map, 0, fout);
  for (uint64_t t = 1; t < cex.size(); ++t ) {
    fout << "#" << t << std::endl;
    dump_diff(cex.at(t), hash_to_value_map, t, fout);
  }
  fout << "#" << cex.size() << std::endl;
}


void VCDWitnessPrinter::DumpTraceToFile(const std::string & vcd_file_name, 
  const std::vector<smt::UnorderedTermMap> & cex) const {
  
  std::ofstream fout(vcd_file_name);
  if (!fout.is_open())
    throw CosaException("Unable to write to : " + vcd_file_name);
  GenHeader(fout);
  DumpValues(fout, cex);
  logger.log(0, "Trace written to " + vcd_file_name);
}

} // namespace cosa
