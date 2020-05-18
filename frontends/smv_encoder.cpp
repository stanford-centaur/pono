#include "smv_encoder.h"
#include "bmc.h"
#include "bmc_simplepath.h"
#include "defaults.h"
#include "interpolantmc.h"
#include "kinduction.h"
#include "prop.h"

using namespace smt;
using namespace cosa;
using namespace std;

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
//parse string
int cosa::SMVEncoder::parseString(std::string newline){
  std::istringstream iss(newline);
  std::istream &s(iss);
  cosa::SMVscanner smvscanner(*this);
  smvscanner.switch_streams(&s);
  cosa::smvparser parse (smvscanner,*this);
  return parse();
}

void cosa::SMVEncoder::processCase(){
    Engine e = BMC;
    unsigned int bound = 10;
    unsigned int verbosity = default_verbosity;
    std::shared_ptr<Prover> prov;
    std::future_status status;
    for (std::pair<int, smt::Term> element : casecheck_)
	  {
        while (transterm_.front().first < element.first){
            rts_.constrain_trans(transterm_.front().second);
            transterm_.pop_front();
        }
        //std::cout << "start case check" <<std::endl; 
        solver_ ->push();
        Property p1(rts_, element.second);
        
        prov = std::make_shared<Bmc>(p1, solver_);
        auto fut = std::async (launch::async,[](std::shared_ptr<Prover> prov,unsigned int bound){return prov->check_until(bound);},prov,bound); 
        std::future_status status;
        status = fut.wait_for(std::chrono::seconds(180));
        while(status != std::future_status::timeout){
          if(status == std::future_status::ready){
            if(fut.get()== UNKNOWN || fut.get()== TRUE){
                 //std::cout<<"find true case" << element.first<<std::endl;
                 rts_.constrain_trans(casestore_[element.first]);
                 break;
            }
            else{
                //std::cout<<"line num" << element.first;
                throw CosaException("case error");
            }  
          }
        }
        if(status == std::future_status::timeout) throw CosaException("case timeout check error"); 
        //std::cout<<"find true case and jump out of the loop" <<std::endl;  
	  }
    solver_->pop();
        while (transterm_.size()>0){
            rts_.constrain_trans(transterm_.front().second);
            transterm_.pop_front();
        }    
}

void cosa::SMVEncoder::preprocess(){
    for (std::pair<int, smt::Term> element : casecheck_)
	  {
        while (transterm_.front().first < element.first){
            rts_.constrain_trans(transterm_.front().second);
            transterm_.pop_front();
        }
        rts_.constrain_trans(casestore_[element.first]);
    }

        while (transterm_.size()>0){
            rts_.constrain_trans(transterm_.front().second);
            transterm_.pop_front();
        }   
}