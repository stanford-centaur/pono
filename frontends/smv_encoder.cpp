#include "smv_encoder.h"

using namespace smt;
using namespace pono;
using namespace std;

int pono::SMVEncoder::parse(std::string filename)
{
  std::ifstream ifs;
  ifs.open(filename);
  if (!ifs.good()) {
    std::cout << "NO input file!" << std::endl;
    exit(-1);
  }
  std::istream & s = ifs;
  pono::SMVscanner smvscanner(*this);
  smvscanner.switch_streams(&s);
  pono::smvparser parse(smvscanner, *this);
  return parse();
}
// parse string
smt::Term pono::SMVEncoder::parseString(std::string newline)
{
  std::istringstream iss(newline);
  std::istream & s(iss);
  pono::SMVscanner smvscanner(*this);
  smvscanner.switch_streams(&s);
  pono::smvparser parse(smvscanner, *this);
  parse();
  return parse_term;
}
// case condition check preprocess
void pono::SMVEncoder::processCase()
{
  std::future_status status;
  solver_->push();
  for (int i = 0; i < casecheck_.size(); i++) {
    Term bad_ = solver_->make_term(smt::PrimOp::Not, casecheck_[i]);
    solver_->assert_formula(bad_);
    auto fut = std::async(
        launch::async,
        [](smt::SmtSolver solver_) {
          Result r = solver_->check_sat();
          return r.is_unsat();
        },
        solver_);
    std::future_status status;
    status = fut.wait_for(std::chrono::seconds(5));
    while (status != std::future_status::timeout) {
      if (status == std::future_status::ready) {
        if (fut.get()) {
          rts_.constrain_trans(casestore_[i]);
          break;
        } else {
          throw PonoException("case error");
        }
      }
    }
    if (status == std::future_status::timeout)
      throw PonoException("case timeout check error");
    solver_->pop();
  }
}
//change the input stream to output stringstream 
int pono::SMVEncoder::parse_flat(std::istream & s)
{
  pono::SMVscanner smvscanner(*this);
  smvscanner.switch_streams(&s);
  pono::smvparser parse(smvscanner, *this);
  return parse();
}

//modular SMV preprocess
std::stringstream pono::SMVEncoder::preprocess()
{
  module_node * main_n;
  if (module_list.find("main") != module_list.end()) {
    main_n = module_list.find("main")->second;
  } else {
    throw PonoException("no main module found");
  }
  std::stringstream str;
  str << "MODULE main" << std::endl;
  main_n->process_main(module_list, str);
  parse_flat(str);
  return str;
}
