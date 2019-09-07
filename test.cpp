#include "frontends/btor2_encoder.h"

#include "smt-switch/boolector_factory.h"

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
  FunctionalTransitionSystem fts(s);

  BTOR2Encoder(filename, fts);

  cout << "Created FunctionalTransitionSystem with:\n";
  cout << "\t" << fts.inputs().size() << " input variables." << endl;
  cout << "\t" << fts.states().size() << " state variables." << endl;

  cout << "Inputs:" << endl;
  for (auto i : fts.inputs())
  {
    cout << "\t" << i << endl;
  }

  cout << "States:" << endl;
  for (auto s : fts.states())
  {
    cout << "\t" << s << endl;
  }

  return 0;
}
