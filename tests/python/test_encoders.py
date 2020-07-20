import pytest
import smt_switch as ss
import pono as c

try:
    import coreir
    COREIR_AVAILABLE=True
except:
    COREIR_AVAILABLE=False

from pathlib import Path


@pytest.mark.skipif(not COREIR_AVAILABLE, reason="Requires coreir python bindings")
@pytest.mark.parametrize("create_solver", ss.solvers.values())
def test_coreir(create_solver):
    coreir_inputs_directory = Path("../encoders/inputs/coreir")
    for f in coreir_inputs_directory.glob("*.json"):
        print("Encoding:", f)
        context = coreir.Context()
        top_mod = context.load_from_file(str(f))

        solver = create_solver(False)
        rts = c.RelationalTransitionSystem(solver)
        ce = c.CoreIREncoder(top_mod, rts)
        t = solver.make_term(True)
        assert ce.trans != t, "Expecting a non-empty transition relation"
