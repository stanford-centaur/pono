#include <unordered_map>
#include <vector>

#include "coreir.h"

#include "core/rts.h"
#include "smt-switch/smt.h"

namespace pono {
class CoreIREncoder
{
 public:
  CoreIREncoder(std::string filename,
                RelationalTransitionSystem & ts,
                bool force_abstract_clock = false)
      : ts_(ts),
        solver_(ts.solver()),
        c_(CoreIR::newContext()),
        num_clocks_(0),
        can_abstract_clock_(true),
        force_abstract_clock_(force_abstract_clock)
  {
    c_->getLibraryManager()->loadLib("commonlib");
    bvsort1_ = solver_->make_sort(smt::BV, 1);
    boolsort_ = solver_->make_sort(smt::BOOL);
    bv1_ = solver_->make_term(1, bvsort1_);
    top_ = read_coreir_file(c_, filename);
    encode();
  }

  CoreIREncoder(CoreIR::Module * m,
                RelationalTransitionSystem & ts,
                bool force_abstract_clock = false)
      : top_(m),
        ts_(ts),
        solver_(ts.solver()),
        c_(m->getContext()),
        num_clocks_(0),
        can_abstract_clock_(true),
        force_abstract_clock_(force_abstract_clock)
  {
    c_->getLibraryManager()->loadLib("commonlib");
    bvsort1_ = solver_->make_sort(smt::BV, 1);
    boolsort_ = solver_->make_sort(smt::BOOL);
    bv1_ = solver_->make_term(1, bvsort1_);

    encode();
  }

 protected:
  static CoreIR::Module * read_coreir_file(CoreIR::Context * c,
                                           std::string filename);
  // encodes the Module * stored in top_
  // to the transition system ts_
  void encode();

  /** looks up the smt encoding for the inst
   *  computes the new term and uses wire_connection
   *  to cache the connections.
   *  Thus if inst.out drives a.in0 and b.in1, it will cache
   *  a.in0 and b.in1 as the term for inst.out in w2term_
   *
   *  @param inst the instance to process
   *  @return the instance output
   */
  CoreIR::Wireable * process_instance(CoreIR::Instance * inst);

  /** computes the next state updates for state elements
   *  this is done as a second pass after all other
   *  instances have been processed so that the drivers
   *  for this state element are already available
   *
   *  relies on can_abstract_clock_ to determine which
   *  encoding to use
   *
   *  @param st the state element Instance *
   */
  void process_state_element(CoreIR::Instance * st);

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

  RelationalTransitionSystem & ts_;
  smt::SmtSolver solver_;
  CoreIR::Context * c_;
  CoreIR::Module * top_;
  size_t num_clocks_;
  bool can_abstract_clock_;  ///< stays true if it's safe to abstract the clock
  bool force_abstract_clock_;  ///< force the clock to be abstracted
                               ///< (synchronizes async behavior)

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
}  // namespace pono
