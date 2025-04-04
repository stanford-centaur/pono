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

#include <functional>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <vector>

#include "core/ts.h"
#include "gmpxx.h"
#include "smt-switch/smt.h"
#include "utils/logger.h"

namespace pono {

struct VCDSignal
{
  std::string vcd_name;  // maybe you want to add this : [N:0]
  std::string full_name;
  std::string hash;
  smt::Term ast;
  uint64_t data_width;
  VCDSignal(const std::string & _vcd_name,
            const std::string & _full_name,
            const std::string & _hash,
            const smt::Term & _ast,
            uint64_t w)
      : vcd_name(_vcd_name),
        full_name(_full_name),
        hash(_hash),
        ast(_ast),
        data_width(w)
  {
  }
};

struct VCDArray : public VCDSignal
{
  std::unordered_map<std::string, std::string> indices2hash;
  VCDArray(const std::string & _vcd_name,
           const std::string & _full_name,
           const smt::Term & _ast,
           uint64_t w)
      : VCDSignal(_vcd_name, _full_name, "", _ast, w)
  {
  }
};

struct VCDScope
{
  std::map<std::string, VCDScope> subscopes;
  std::map<std::string, VCDSignal> wires;
  std::map<std::string, VCDSignal> regs;
  std::map<std::string, VCDArray> arrays;
};  // struct VCDScope

class VCDWitnessPrinter
{
 public:
  // types
  typedef std::map<std::string, std::vector<std::string>> per_mem_indices;

 protected:
  const smt::UnorderedTermSet & inputs_;
  const smt::UnorderedTermSet & states_;
  const std::unordered_map<std::string, smt::Term> & named_terms_;
  const std::vector<smt::UnorderedTermMap> & cex_;

  VCDScope root_scope_;
  // hash id --> signal object
  std::vector<VCDSignal *> allsig_bv_;
  std::vector<VCDArray *> allsig_array_;

  // given a name like a.b.c, find the right scope and
  // create if it does not exists
  void check_insert_scope(std::string full_name,
                          bool is_reg,
                          const smt::Term & ast);
  // another function for array maybe?
  void check_insert_scope_array(std::string full_name,
                                const std::unordered_set<std::string> & indices,
                                bool has_default,
                                const smt::Term & ast);

  uint64_t hash_id_cnt_;
  std::string new_hash_id();

  uint64_t property_id_cnt_;
  std::string new_property_id();

  void dump_current_scope(std::ostream & fout, const VCDScope *) const;

  void dump_all(const smt::UnorderedTermMap & valmap,
                std::unordered_map<std::string, std::string> & valbuf,
                uint64_t tick,
                std::ostream & fout) const;
  void dump_diff(const smt::UnorderedTermMap & valmap,
                 std::unordered_map<std::string, std::string> & valprev,
                 uint64_t tick,
                 std::ostream & fout) const;

 protected:
  void DumpScopes(std::ostream & fout) const;
  void GenHeader(std::ostream & fout) const;
  void DumpValues(std::ostream & fout) const;

 public:
  VCDWitnessPrinter(
      const TransitionSystem & ts,
      const std::vector<smt::UnorderedTermMap> & cex,
      const std::unordered_map<std::string, std::string> & symbol_map = {});

  void dump_trace_to_file(const std::string & vcd_file_name) const;
  void debug_dump() const;

};  // class VCDWitnessPrinter

}  // namespace pono
