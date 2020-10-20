import os
import pytest
import smt_switch as ss
import pono

try:
    import coreir
    COREIR_AVAILABLE=True
except:
    COREIR_AVAILABLE=False

from pathlib import Path


@pytest.mark.skipif(not COREIR_AVAILABLE, reason="Requires coreir python bindings")
@pytest.mark.skipif(not hasattr(pono, "CoreIREncoder"), reason="Requires building with coreir")
@pytest.mark.parametrize("create_solver", ss.solvers.values())
def test_coreir(create_solver):
    file_dir = Path(os.path.dirname(__file__))
    coreir_inputs_directory = file_dir / Path("../encoders/inputs/coreir")
    coreir_files = list(coreir_inputs_directory.glob("*.json"))
    assert coreir_files
    for f in coreir_files:
        print("Encoding:", f)
        context = coreir.Context()
        context.load_library("commonlib")
        top_mod = context.load_from_file(str(f))

        solver = create_solver(False)
        rts = pono.RelationalTransitionSystem(solver)
        ce = pono.CoreIREncoder(top_mod, rts)
        t = solver.make_term(True)
        assert rts.trans != t, "Expecting a non-empty transition relation"
