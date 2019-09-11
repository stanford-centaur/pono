#include "frontends/btor2_encoder.h"

#include "smt-switch/boolector_factory.h"

#include "prop.h"

#include "bmc.h"

#include <iostream>

using namespace cosa;
using namespace smt;
using namespace std;

int main(int argc, char ** argv)
{
  if (argc != 2)
  {
    throw "Expecting btor file as input.";
  }

  string filename = argv[1];

  SmtSolver s = BoolectorSolverFactory::create();
  s->set_opt("produce-models", true);
  RelationalTransitionSystem rts(s);

  BTOR2Encoder btor_enc(filename, rts);

  cout << "Created FunctionalTransitionSystem with:\n";
  cout << "\t" << rts.inputs().size() << " input variables." << endl;
  cout << "\t" << rts.states().size() << " state variables." << endl;

  cout << "Inputs:" << endl;
  for (auto i : rts.inputs())
  {
    cout << "\t" << i << endl;
  }

  cout << "States:" << endl;
  for (auto s : rts.states())
  {
    cout << "\t" << s << endl;
  }

  cout << "Init:" << endl;
  cout << "\t" << rts.init() << endl;

  cout << "Trans:" << endl;
  cout << "\t" << rts.trans() << endl;

  for (auto i : btor_enc.badvec()) {
    Property p(rts, s->make_term(PrimOp::Not, i));
    cout << "Property:" << endl;
    cout << p.prop() << endl;

    Bmc bmc(p, s);

    ProverResult r = bmc.check_until(20);
    if (r == FALSE) {
      cout << i << " is FALSE" << endl;
    } else if (r == TRUE) {
      cout << i << " is TRUE" << endl;
    } else {
      cout << i << " is UNKNOWN" << endl;
    }
  }
  
  return 0;
}
