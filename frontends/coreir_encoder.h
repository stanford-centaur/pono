#include <unordered_map>
#include <vector>

#include "coreir.h"

#include "core/fts.h"
#include "smt-switch/smt.h"

namespace cosa {
class CoreIREncoder
{
 public:
  CoreIREncoder(std::string filename, TransitionSystem & ts)
      : ts_(ts),
        solver_(ts.solver()),
        c_(CoreIR::newContext()),
        num_clocks_(0),
        can_abstract_clock_(true)
  {
    parse(filename);
  }

 protected:
  static CoreIR::Module * read_coreir_file(CoreIR::Context * c,
                                           std::string filename);
  void parse(std::string filename);
  void process_instance(CoreIR::Instance * inst);

  TransitionSystem & ts_;
  smt::SmtSolver solver_;
  CoreIR::Context * c_;
  CoreIR::Module * top_;
  CoreIR::Module * mod_;
  CoreIR::ModuleDef * def_;
  size_t num_clocks_;
  bool can_abstract_clock_;

  // conversion data structures
  std::unordered_map<CoreIR::Wireable *, smt::Term> w2term_;

  // useful temporary variables
  smt::Term t_;
  smt::Sort sort_;
  CoreIR::Instance * inst_;
  CoreIR::Type * type_;

  const std::vector<std::string> passes_ = {
    "clockifyinterface",
    "rungenerators",
    "removeconstduplicates",
    "deletedeadinstances",
    "rungenerators",
    "cullgraph",
    "removebulkconnections",
    "removeunconnected",
    "flatten",
    // might not be necessary if it's already flat
    "flattentypes",
    "packconnections",
    "cullzexts"
  };
};
}  // namespace cosa
