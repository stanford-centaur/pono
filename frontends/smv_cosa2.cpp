#include "node.h"
#include "rts.h"
#include "smt-switch/boolector_factory.h"
#include "smt-switch/smt.h"
#include "smv_encoder.h"

#include "bmc.h"
#include "bmc_simplepath.h"
#include "defaults.h"
#include "interpolant.h"
#include "kinduction.h"
#include "prop.h"
#include "logger.h"
//#include "smt-switch/msat_factory.h"

using namespace smt;
using namespace cosa;
using namespace std;
//


int main(int, char** argv) {
  
  SmtSolver s = BoolectorSolverFactory::create();
  //SmtSolver s = MsatSolverFactory::create();
        s->set_opt("produce-models", "true");
      s->set_opt("incremental", "true");
  RelationalTransitionSystem rts(s);
  std::cout << "++++++++++++++++++++++++++ PRINT SIMPLE SYSTEM ++++++++++++++++++++++++++" << std::endl;
  std::cout << argv[1] << std::endl;
  smvEncoder smv_enc(argv[1], rts);
  unsigned int prop_idx = default_prop_idx;
  //std::cout << "++++++++++++++++++++++++++ PRINT SIMPLE SYSTEM ++++++++++++++++++++++++++" << std::endl;
  std::cout << "INPUTS" << std::endl;
  for (auto i : rts.inputs())
  {
    std::cout << "\t" << i << std::endl;
  }
  std::cout << "STATES" << std::endl;
  for (auto s : rts.states())
  {
    std::cout << "\t" << s << std::endl;
  }

  std::cout << "INIT " << rts.init() << std::endl;
  std::cout << "TRANS " << rts.trans() << std::endl;
  //std::cout << "INVARSPEC " << smv_enc.propvec_[0] << std::endl;

    Engine engine = BMC;
    unsigned int bound = 10;
    unsigned int verbosity = default_verbosity;
    Property p(rts, smv_enc.propvec_[0]);
    logger.log(1, "Solving property: {}", p.prop());
  std::cout << "INVARSPEC " << p.prop() << std::endl;
    std::shared_ptr<Prover> prover;
    if (engine == BMC) {
      prover = std::make_shared<Bmc>(p, s);
    } else if (engine == BMC_SP) {
      prover = std::make_shared<BmcSimplePath>(p, s);
    } else if (engine == KIND) {
      prover = std::make_shared<KInduction>(p, s);
    } else {
      throw CosaException("Unimplemented engine.");
    }

     ProverResult r = prover->check_until(bound);
    if (r == FALSE) {
      cout << "Find a counterexample" << endl;
      cout << "b" << prop_idx << endl;
      vector<UnorderedTermMap> cex;
      if (prover->witness(cex)) {
        //print_witness_btor(smv_enc, cex);
      }
      return 1;
    } else if (r == TRUE) {
      //cout << "unsat" << endl;
      cout<<"Proved"<<endl;
      cout << "b" << prop_idx << endl;
      return 0;
    } else {
      cout << "unknown" << endl;
      cout << "b" << prop_idx << endl;
      return 2;
    }

  return 3;
}
