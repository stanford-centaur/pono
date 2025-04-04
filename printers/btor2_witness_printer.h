/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Makai Mann, Ahmed Irfan
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

#include <iterator>
#include <map>
#include <sstream>
#include <vector>

#include "gmpxx.h"
#include "smt-switch/smt.h"
#include "utils/logger.h"
#include "utils/str_util.h"

namespace pono {

std::string as_bits(std::string val)
{
  // TODO: this makes assumptions on format of value from boolector
  //       to support other solvers, we need to be more general
  std::string res = val;

  if (val.length() < 2) {
    throw PonoException("Don't know how to interpret value: " + val);
  }

  if (res.substr(0, 2) == "#b") {
    // remove the #b prefix
    res = res.substr(2, val.length() - 2);
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

// Returns true iff 'term' appears in either the state variables or in
// the input variables of 'ts'. If COI is applied, then 'ts' might
// have different state and input variables than the original
// transition system that was parsed by the BTOR encoder.
// NOTE: we call this check also when COI is not applied, hence some
// overhead will occur in witness printing.
bool appears_in_ts_coi(const smt::Term & term, const TransitionSystem & ts)
{
  const auto & it_states = ts.statevars().find(term);
  if (it_states != ts.statevars().end()) return true;

  const auto & it_inputs = ts.inputvars().find(term);
  if (it_inputs != ts.inputvars().end()) return true;

  return false;
}

void print_btor_vals_at_time(
    const smt::TermVec & vec,
    const smt::UnorderedTermMap & valmap,
    unsigned int time,
    const TransitionSystem & ts,
    const std::unordered_map<std::string, std::string> & symbol_map)
{
  smt::SortKind sk;
  smt::TermVec store_children(3);
  for (size_t i = 0, size = vec.size(); i < size; ++i) {
    // Do not print if term 'vec[i]' does not appear in COI. When not
    // using COI, this check always returns true.
    if (!appears_in_ts_coi(vec[i], ts)) continue;
    sk = vec[i]->get_sort()->get_sort_kind();
    if (sk == smt::BV) {
      // TODO: this makes assumptions on format of value from boolector
      //       to support other solvers, we need to be more general
      logger.log(0,
                 "{} {} {}@{}",
                 i,
                 as_bits(valmap.at(vec[i])->to_string()),
                 lookup_or_key(symbol_map, vec[i]->to_string()),
                 time);
    } else if (sk == smt::ARRAY) {
      smt::Term tmp = valmap.at(vec[i]);
      while (tmp->get_op() == smt::Store) {
        int num = 0;
        for (auto c : tmp) {
          store_children[num] = c;
          num++;
        }

        logger.log(0,
                   "{} [{}] {} {}@{}",
                   i,
                   as_bits(store_children[1]->to_string()),
                   as_bits(store_children[2]->to_string()),
                   lookup_or_key(symbol_map, vec[i]->to_string()),
                   time);
        tmp = store_children[0];
      }

      if (tmp->get_op().is_null()
          && tmp->get_sort()->get_sort_kind() == smt::ARRAY) {
        smt::Term const_val = *(tmp->begin());
        logger.log(0,
                   "{} {} {}@{}",
                   i,
                   as_bits(const_val->to_string()),
                   lookup_or_key(symbol_map, vec[i]->to_string()),
                   time);
      }

    } else {
      throw PonoException("Unhandled sort kind: " + ::smt::to_string(sk));
    }
  }
}

void print_btor_vals_at_time(
    const std::map<uint64_t, smt::Term> m,
    const smt::UnorderedTermMap & valmap,
    unsigned int time,
    const TransitionSystem & ts,
    const std::unordered_map<std::string, std::string> & symbol_map)
{
  smt::SortKind sk;
  smt::TermVec store_children(3);
  for (auto entry : m) {
    // Do not print if term 'entry.second' does not appear in COI. When not
    // using COI, this check always returns true.
    if (!appears_in_ts_coi(entry.second, ts)) continue;
    sk = entry.second->get_sort()->get_sort_kind();
    if (sk == smt::BV) {
      // TODO: this makes assumptions on format of value from boolector
      //       to support other solvers, we need to be more general
      // Remove the #b prefix
      logger.log(0,
                 "{} {} {}@{}",
                 entry.first,
                 as_bits(valmap.at(entry.second)->to_string()),
                 lookup_or_key(symbol_map, entry.second->to_string()),
                 time);
    } else if (sk == smt::ARRAY) {
      smt::Term tmp = valmap.at(entry.second);
      while (tmp->get_op() == smt::Store) {
        int num = 0;
        for (auto c : tmp) {
          store_children[num] = c;
          num++;
        }

        logger.log(0,
                   "{} [{}] {} {}@{}",
                   entry.first,
                   as_bits(store_children[1]->to_string()),
                   as_bits(store_children[2]->to_string()),
                   lookup_or_key(symbol_map, entry.second->to_string()),
                   time);
        tmp = store_children[0];
      }

      if (tmp->get_op().is_null()
          && tmp->get_sort()->get_sort_kind() == smt::ARRAY) {
        smt::Term const_val = *(tmp->begin());
        logger.log(0,
                   "{} {} {}@{}",
                   entry.first,
                   as_bits(const_val->to_string()),
                   lookup_or_key(symbol_map, entry.second->to_string()),
                   time);
      }

    } else {
      throw PonoException("Unhandled sort kind: " + ::smt::to_string(sk));
    }
  }
}

void print_witness_btor(const BTOR2Encoder & btor_enc,
                        const std::vector<smt::UnorderedTermMap> & cex,
                        const TransitionSystem & ts)
{
  const smt::TermVec inputs = btor_enc.inputsvec();
  const smt::TermVec states = btor_enc.statesvec();
  const std::map<uint64_t, smt::Term> no_next_states =
      btor_enc.no_next_statevars();
  bool has_states_without_next = no_next_states.size();

  logger.log(0, "#0");
  print_btor_vals_at_time(states, cex.at(0), 0, ts, btor_enc.get_symbol_map());

  for (size_t k = 0, cex_size = cex.size(); k < cex_size; ++k) {
    // states without next
    if (k && has_states_without_next) {
      logger.log(0, "#{}", k);
      print_btor_vals_at_time(
          no_next_states, cex.at(k), k, ts, btor_enc.get_symbol_map());
    }

    // inputs
    if (k < cex_size) {
      logger.log(0, "@{}", k);
      print_btor_vals_at_time(
          inputs, cex.at(k), k, ts, btor_enc.get_symbol_map());
    }
  }

  logger.log(0, ".");
}

}  // namespace pono
