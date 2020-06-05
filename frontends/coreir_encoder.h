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
    bvsort1_ = solver_->make_sort(smt::BV, 1);
    boolsort_ = solver_->make_sort(smt::BOOL);
    bv1_ = solver_->make_term(1, bvsort1_);
    parse(filename);
  }

 protected:
  static CoreIR::Module * read_coreir_file(CoreIR::Context * c,
                                           std::string filename);
  void parse(std::string filename);
  CoreIR::Wireable * process_instance(CoreIR::Instance * inst);
  /** wires up a connection by modifying w2term_
   *  to point to the correct term
   *  if necessary, it will create a "forward reference"
   *  for the destination, e.g. if the destination is
   *  connected via a select (and not the whole wire)
   *  @param conn the connection
   */
  void wire_connection(CoreIR::Connection conn);
  /** Gets sort for an arbitrary Wireable using the Type
   *  @param w the Wireable
   *  @return an Smt-Switch Bool or BV sort
   *  Note: currently does not support memories
   */
  smt::Sort compute_sort(CoreIR::Wireable * w);

  TransitionSystem & ts_;
  smt::SmtSolver solver_;
  CoreIR::Context * c_;
  CoreIR::Module * top_;
  size_t num_clocks_;
  bool can_abstract_clock_;

  // conversion data structures
  std::unordered_map<CoreIR::Wireable *, smt::Term> w2term_;

  // useful reusable variables
  smt::Sort bvsort1_;
  smt::Sort boolsort_;
  smt::Term bv1_;

  // useful temporary variables
  smt::Term t_;
  smt::Sort sort_;
  CoreIR::Instance * inst_;
  CoreIR::Type * type_;
  CoreIR::Module * mod_;
  CoreIR::ModuleDef * def_;

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
