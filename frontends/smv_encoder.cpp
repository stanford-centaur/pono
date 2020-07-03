#include "smv_encoder.h"

using namespace smt;
using namespace pono;
using namespace std;

int pono::SMVEncoder::parse(std::string filename)
{
  std::ifstream ifs;
  ifs.open(filename);
  // make sure file is valid:
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
int pono::SMVEncoder::parseString(std::string newline)
{
  std::istringstream iss(newline);
  std::istream & s(iss);
  pono::SMVscanner smvscanner(*this);
  smvscanner.switch_streams(&s);
  pono::smvparser parse(smvscanner, *this);
  return parse();
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
    if (status == std::future_status::timeout) {
      throw PonoException("case timeout check error");
    }
  }
  solver_->pop();
}