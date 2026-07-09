from __future__ import annotations

from typing import Callable

import pono
import pytest
import smt_switch as ss


@pytest.mark.parametrize("create_solver", ss.solvers.values())
def test_replace_terms(create_solver: Callable[[bool], ss.SmtSolver]) -> None:
    solver = create_solver(create_solver is ss.solvers.get("yices2"))
    bvsort8 = solver.make_sort(ss.sortkinds.BV, 8)
    fts = pono.FunctionalTransitionSystem(solver)
    x = fts.make_statevar("x", bvsort8)
    a = fts.make_statevar("a", bvsort8)
    b = fts.make_statevar("b", bvsort8)

    amb = fts.make_term(ss.primops.BVMul, a, b)
    fts.assign_next(x, amb)

    freshinput = fts.make_inputvar("freshinput", bvsort8)

    to_replace = {amb: freshinput}
    fts.replace_terms(to_replace)

    assert fts.state_updates[x] == freshinput, "Expecting terms to be replaced"
