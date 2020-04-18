import pytest
import smt_switch as ss
import cosa2 as c

def build_simple_ts(solver, TS):
    bvsort = solver.make_sort(ss.sortkinds.BV, 8)

    ts = TS(solver)
    x = ts.make_state('x', bvsort)
    y = ts.make_state('y', bvsort)
    xp1 = solver.make_term(ss.primops.BVAdd, x, solver.make_term(1, bvsort))
    ts.name_term('xp1', xp1)
    ts.assign_next(x, xp1)
    assert ts.state_updates[x] == xp1

    return solver, ts


@pytest.mark.parametrize("create_solver", ss.solvers.values())
def test_cons_fts(create_solver):
    solver = create_solver()
    solver, ts = build_simple_ts(solver, c.FunctionalTransitionSystem)

@pytest.mark.parametrize("create_solver", ss.solvers.values())
def test_query_fts(create_solver):
    solver = create_solver()
    solver, ts = build_simple_ts(solver, c.FunctionalTransitionSystem)

    assert len(ts.states) == 2
    assert len(ts.state_updates) == 1
    assert len(ts.named_terms) == 1
    assert ts.is_functional()

    states = list(ts.states)
    try:
        ts.constrain_trans(solver.make_term(ss.primops.Equal, ts.next(states[0]), ts.next(states[1])))
        assert False
    except Exception as e:
        pass

@pytest.mark.parametrize("create_solver", ss.solvers.values())
def test_func_update_fts(create_solver):
    solver = create_solver()
    solver, ts = build_simple_ts(solver, c.FunctionalTransitionSystem)

    states = list(ts.states)
    try:
        ts.assign_next(states[0], ts.next(states[1]))
        assert False
    except:
        pass

@pytest.mark.parametrize("create_solver", ss.solvers.values())
def test_cons_rts(create_solver):
    solver = create_solver()
    solver, ts = build_simple_ts(solver, c.RelationalTransitionSystem)

@pytest.mark.parametrize("create_solver", ss.solvers.values())
def test_query_rts(create_solver):
    solver = create_solver()
    solver, ts = build_simple_ts(solver, c.RelationalTransitionSystem)

    assert len(ts.states) == 2
    assert len(ts.state_updates) == 1
    assert len(ts.named_terms) == 1
    assert not ts.is_functional()

    states = list(ts.states)
    try:
        ts.constrain_trans(solver.make_term(ss.primops.Equal, ts.next(states[0]), ts.next(states[1])))
    except Exception as e:
        assert False


