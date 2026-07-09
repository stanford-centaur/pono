from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Callable

import pono
import pytest
import smt_switch as ss


@pytest.mark.parametrize("create_solver", ss.solvers.values())
def test_vcd_trace(create_solver: Callable[[bool], ss.SmtSolver]) -> None:
    if create_solver is ss.solvers.get("z3"):
        pytest.skip(reason="Z3's hexadecimal format is not supported yet")
    solver = create_solver(create_solver is ss.solvers.get("yices2"))
    solver.set_opt("incremental", "true")
    solver.set_opt("produce-models", "true")
    if "btor" in ss.solvers and ss.solvers["btor"] is create_solver:
        solver.set_opt("base-context-1", "true")
    bvsort8 = solver.make_sort(ss.sortkinds.BV, 8)
    ts = pono.FunctionalTransitionSystem(solver)
    x = ts.make_statevar("x", bvsort8)
    ts.constrain_init(
        solver.make_term(ss.primops.Equal, x, solver.make_term(0, x.get_sort()))
    )
    ts.assign_next(
        x, ts.make_term(ss.primops.BVAdd, x, solver.make_term(1, x.get_sort()))
    )

    prop_term = solver.make_term(ss.primops.BVUle, x, solver.make_term(9, x.get_sort()))

    prop = pono.Property(solver, prop_term)
    bmc = pono.Bmc(prop, ts, solver)
    res = bmc.check_until(10)
    assert not res, "res should be false, not just unknown (i.e. None)"

    witness = bmc.witness()

    with tempfile.NamedTemporaryFile() as temp:
        temp_path = Path(temp.name)
        assert temp_path.stat().st_size == 0, "Expect file to start empty"
        vcd_printer = pono.VCDWitnessPrinter(ts, witness)
        vcd_printer.dump_trace_to_file(temp.name)
        assert temp_path.stat().st_size, "Expect file to be non-empty"
