/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Makai Mann, Ahmed Irfan, Po-Chun Chien, Áron Ricardo Perez-Lopez
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
#include "printers/btor2_witness_printer.h"

#include <cstdint>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "core/ts.h"
#include "frontends/btor2_encoder.h"
#include "gmpxx.h"
#include "smt-switch/smt.h"
#include "utils/exceptions.h"
#include "utils/logger.h"
#include "utils/str_util.h"

namespace pono {

std::string as_bits(const std::uint64_t value, const std::uint64_t bit_width)
{
  assert(bit_width > 0);
  std::string res(bit_width, '0');
  res.reserve(bit_width);

  std::uint64_t bit_mask = 1;
  for (std::uint64_t i = 0; i < bit_width; i++) {
    if (value & bit_mask) {
      res[bit_width - 1 - i] = '1';
    }
    bit_mask <<= 1;
  }
  return res;
}

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
    std::size_t width = std::stoull(width_str);
    mpz_class cval(res);
    res = cval.get_str(2);
    std::size_t strlen = res.length();

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

void print_btor_val_at_time(
    std::size_t btor_id,
    smt::Term term,
    unsigned int time,
    const TransitionSystem & ts,
    const smt::UnorderedTermMap & valmap,
    const std::unordered_map<std::string, std::string> & symbol_map,
    std::ostream & output_stream)
{
  // Do not print if term term does not appear in COI. When not
  // using COI, this check always returns true.
  if (!appears_in_ts_coi(term, ts)) {
    return;
  }
  smt::SortKind sk = term->get_sort()->get_sort_kind();
  if (sk == smt::BV) {
    // TODO: this makes assumptions on format of value from boolector
    //       to support other solvers, we need to be more general
    logger.log_to_stream(0,
                         output_stream,
                         "{} {} {}@{}",
                         btor_id,
                         as_bits(valmap.at(term)->to_string()),
                         lookup_or_key(symbol_map, term->to_string()),
                         time);
  } else if (sk == smt::ARRAY) {
    smt::Term tmp = valmap.at(term);
    smt::TermVec store_children(3);
    std::unordered_set<std::uint64_t> written_ids;
    while (tmp->get_op() == smt::Store) {
      int num = 0;
      for (auto c : tmp) {
        store_children[num] = c;
        num++;
      }
      logger.log_to_stream(0,
                           output_stream,
                           "{} [{}] {} {}@{}",
                           btor_id,
                           as_bits(store_children[1]->to_string()),
                           as_bits(store_children[2]->to_string()),
                           lookup_or_key(symbol_map, term->to_string()),
                           time);
      written_ids.insert(store_children[1]->to_int());
      tmp = store_children[0];
    }
    if (tmp->get_op().is_null()
        && tmp->get_sort()->get_sort_kind() == smt::ARRAY) {
      smt::Term const_val = *(tmp->begin());
      if (const_val->to_int() == 0) {
        // Skip when value is 0.
        // When unspecified, btorsim assumes 0 as the default value.
        // Although the correct semantics treat unspecified values as "don't
        // care", printing 0s for arrays with large index widths would explode
        // the witness size.
        return;
      }
      std::uint64_t index_width =
          term->get_sort()->get_indexsort()->get_width();
      std::uint64_t index_size = 1ULL << index_width;
      for (std::size_t j = 0; j < index_size; ++j) {
        if (written_ids.find(j) != written_ids.end()) {
          continue;
        }
        logger.log_to_stream(0,
                             output_stream,
                             "{} [{}] {} {}@{}",
                             btor_id,
                             as_bits(j, index_width),
                             as_bits(const_val->to_string()),
                             lookup_or_key(symbol_map, term->to_string()),
                             time);
      }
    }
  } else {
    throw PonoException("Unhandled sort kind: " + ::smt::to_string(sk));
  }
}

void print_btor_vals_at_time(
    const smt::TermVec & vec,
    const smt::UnorderedTermMap & valmap,
    unsigned int time,
    const TransitionSystem & ts,
    const std::unordered_map<std::string, std::string> & symbol_map,
    std::ostream & output_stream,
    const smt::UnorderedTermSet & skip_terms = {})
{
  for (std::size_t i = 0, size = vec.size(); i < size; ++i) {
    // Skip terms that are present in the given set
    // (e.g., state variables initialized at #0)
    if (skip_terms.find(vec.at(i)) != skip_terms.end()) {
      continue;
    }
    print_btor_val_at_time(
        i, vec[i], time, ts, valmap, symbol_map, output_stream);
  }
}

void print_btor_vals_at_time(
    const std::map<uint64_t, smt::Term> m,
    const smt::UnorderedTermMap & valmap,
    unsigned int time,
    const TransitionSystem & ts,
    const std::unordered_map<std::string, std::string> & symbol_map,
    std::ostream & output_stream)
{
  for (auto entry : m) {
    print_btor_val_at_time(
        entry.first, entry.second, time, ts, valmap, symbol_map, output_stream);
  }
}

void print_witness_btor(const BTOR2Encoder & btor_enc,
                        const std::vector<smt::UnorderedTermMap> & cex,
                        const TransitionSystem & ts,
                        std::ostream & output_stream)
{
  const smt::TermVec inputs = btor_enc.inputsvec();
  const smt::TermVec states = btor_enc.statesvec();
  const std::map<uint64_t, smt::Term> no_next_states =
      btor_enc.no_next_statevars();
  bool has_states_without_next = !no_next_states.empty();

  logger.log_to_stream(0, output_stream, "#0");
  print_btor_vals_at_time(states,
                          cex.at(0),
                          0,
                          ts,
                          btor_enc.get_symbol_map(),
                          output_stream,
                          btor_enc.initialized_statevars());

  for (std::size_t k = 0, cex_size = cex.size(); k < cex_size; ++k) {
    // states without next
    if (k && has_states_without_next) {
      logger.log_to_stream(0, output_stream, "#{}", k);
      print_btor_vals_at_time(no_next_states,
                              cex.at(k),
                              k,
                              ts,
                              btor_enc.get_symbol_map(),
                              output_stream);
    }

    // inputs
    if (k < cex_size) {
      logger.log_to_stream(0, output_stream, "@{}", k);
      print_btor_vals_at_time(
          inputs, cex.at(k), k, ts, btor_enc.get_symbol_map(), output_stream);
    }
  }

  logger.log_to_stream(0, output_stream, ".");
}

void dump_witness_btor(const BTOR2Encoder & btor_enc,
                       std::vector<smt::UnorderedTermMap> & cex,
                       const TransitionSystem & ts,
                       const unsigned int prop_idx,
                       const std::string & witness_filename)
{
  std::ofstream fout(witness_filename);
  if (!fout.is_open()) {
    throw PonoException("Failed to open file: " + witness_filename);
  }

  fout << "sat" << std::endl;
  fout << "b" << prop_idx << std::endl;
  print_witness_btor(btor_enc, cex, ts, fout);
}

}  // namespace pono
