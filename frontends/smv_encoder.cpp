#include "smv_encoder.h"

int cosa::SMVEncoder::parse(std::string filename){
    std::ifstream ifs;
    ifs.open(filename);
  // make sure file is valid:
    if (!ifs.good()) {
    std::cout << "NO input file!" << std::endl;
    exit(-1);
  }
  std::istream& s= ifs;
  cosa::SMVscanner smvscanner(*this);
  smvscanner.switch_streams(&s);
  cosa::smvparser parse (smvscanner,*this);
  return parse();
}

int cosa::SMVEncoder::parseString(std::string newline){
  std::istringstream iss(newline);
  std::istream &s(iss);
  cosa::SMVscanner smvscanner(*this);
  smvscanner.switch_streams(&s);
  cosa::smvparser parse (smvscanner,*this);
  return parse();
}

void cosa::SMVEncoder::processCase(){
    Engine engine = BMC;
    unsigned int bound = 10;
    unsigned int verbosity = default_verbosity;
    Property p(rts, smc);
    logger.log(1, "Solving property: {}", p.prop());
    std::cout << "INVARSPEC " << p.prop() << std::endl;
    std::shared_ptr<Prover> prover;
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
}