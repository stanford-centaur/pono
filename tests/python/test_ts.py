import pytest
import smt_switch as ss
import pono

def build_simple_ts(solver, TS):
    bvsort = solver.make_sort(ss.sortkinds.BV, 8)

    ts = TS(solver)
    x = ts.make_statevar('x', bvsort)
    y = ts.make_statevar('y', bvsort)
    xp1 = solver.make_term(ss.primops.BVAdd, x, solver.make_term(1, bvsort))
    ts.name_term('xp1', xp1)
    ts.assign_next(x, xp1)
    assert ts.state_updates[x] == xp1

    ts.add_invar(solver.make_term(ss.primops.Equal, y, solver.make_term(4, bvsort)))

    return solver, ts


@pytest.mark.parametrize("create_solver", ss.solvers.values())
def test_cons_fts(create_solver):
    solver = create_solver(False)
    solver, ts = build_simple_ts(solver, pono.FunctionalTransitionSystem)

@pytest.mark.parametrize("create_solver", ss.solvers.values())
def test_query_fts(create_solver):
    solver = create_solver(False)
    solver, ts = build_simple_ts(solver, pono.FunctionalTransitionSystem)

    assert len(ts.statevars) == 2
    assert len(ts.state_updates) == 1
    assert len(ts.named_terms) == 5, "expecting a named term for each curr/next state var and explicitly named term"
    assert len(ts.constraints) == 1, "expecting the added constraint over current state vars only"
    assert ts.is_functional()
    assert not ts.is_deterministic(), "not deterministic because no update for y"

    states = list(ts.statevars)
    try:
        ts.constrain_trans(solver.make_term(ss.primops.Equal, ts.next(states[0]), ts.next(states[1])))
        assert False
    except Exception as e:
        pass

@pytest.mark.parametrize("create_solver", ss.solvers.values())
def test_func_update_fts(create_solver):
    solver = create_solver(False)
    solver, ts = build_simple_ts(solver, pono.FunctionalTransitionSystem)

    states = list(ts.statevars)
    try:
        ts.assign_next(states[0], ts.next(states[1]))
        assert False
    except:
        pass

@pytest.mark.parametrize("create_solver", ss.solvers.values())
def test_cons_rts(create_solver):
    solver = create_solver(False)
    solver, ts = build_simple_ts(solver, pono.RelationalTransitionSystem)

@pytest.mark.parametrize("create_solver", ss.solvers.values())
def test_query_rts(create_solver):
    solver = create_solver(False)
    solver, ts = build_simple_ts(solver, pono.RelationalTransitionSystem)

    assert len(ts.statevars) == 2
    assert len(ts.state_updates) == 1
    assert len(ts.named_terms) == 5, "expecting a named term for each curr/next state var and explicitly named term"
    assert len(ts.constraints) == 1, "expecting the added constraint over current state vars only"
    assert not ts.is_functional()

    states = list(ts.statevars)
    try:
        ts.constrain_trans(solver.make_term(ss.primops.Equal, ts.next(states[0]), ts.next(states[1])))
    except Exception as e:
        assert False

