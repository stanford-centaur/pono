#include "smv_encoder.h"

using namespace smt;
using namespace cosa;
using namespace std;

int cosa::SMVEncoder::parse(std::string filename)
{
  std::ifstream ifs;
  ifs.open(filename);
  if (!ifs.good()) {
    std::cout << "NO input file!" << std::endl;
    exit(-1);
  }
  std::istream & s = ifs;
  cosa::SMVscanner smvscanner(*this);
  smvscanner.switch_streams(&s);
  cosa::smvparser parse(smvscanner, *this);
  return parse();
}
// parse string
int cosa::SMVEncoder::parseString(std::string newline)
{
  std::istringstream iss(newline);
  std::istream & s(iss);
  cosa::SMVscanner smvscanner(*this);
  smvscanner.switch_streams(&s);
  cosa::smvparser parse(smvscanner, *this);
  return parse();
}
// case condition check preprocess
void cosa::SMVEncoder::processCase()
{
  std::future_status status;
  for (std::pair<int, smt::Term> element : casecheck_) {
    solver_->push();
    Term bad_ = solver_->make_term(smt::PrimOp::Not, element.second);
    solver_->assert_formula(bad_);
    auto fut = std::async(
        launch::async,
        [](smt::SmtSolver solver_) {
          Result r = solver_->check_sat();
          return r.is_sat();
        },
        solver_);
    std::future_status status;
    status = fut.wait_for(std::chrono::seconds(5));
    while (status != std::future_status::timeout) {
      if (status == std::future_status::ready) {
        if (fut.get() == false) {
          rts_.constrain_trans(casestore_[element.first]);
          break;
        } else {
          throw CosaException("case error");
        }
      }
    }
    if (status == std::future_status::timeout)
      throw CosaException("case timeout check error");
    solver_->pop();
  }
}
//change the input stream to output stringstream 
int cosa::SMVEncoder::parse_flat(std::istream & s)
{
  cosa::SMVscanner smvscanner(*this);
  smvscanner.switch_streams(&s);
  cosa::smvparser parse(smvscanner, *this);
  return parse();
}

//modular SMV preprocess
void cosa::SMVEncoder::preprocess()
{
  module_node * main_n;
  if (module_list.find("main") != module_list.end()) {
    main_n = module_list.find("main")->second;
  } else {
    throw CosaException("no main module find");
  }
  std::stringstream str;
  str << "MODULE main" << std::endl;
  main_n->pre(module_list, str);
  parse_flat(str);
}