import pytest
import smt_switch as ss
import pycosa2 as c
from typing import Set

def get_free_vars(t:ss.Term)->Set[ss.Term]:
    to_visit = [t]
    visited = set()

    free_vars = set()

    while to_visit:
        t = to_visit[-1]
        to_visit = to_visit[:-1]

        if t in visited:
            continue
        else:
            for tt in t:
                to_visit.append(tt)

            if t.is_symbolic_const():
                free_vars.add(t)

    return free_vars


@pytest.mark.parametrize("create_solver", ss.solvers.values())
def test_fts_unroller(create_solver):
    solver = create_solver()
    bvsort4 = solver.make_sort(ss.sortkinds.BV, 4)
    bvsort8 = solver.make_sort(ss.sortkinds.BV, 8)
    arrsort = solver.make_sort(ss.sortkinds.ARRAY, bvsort4, bvsort8)

    ts = c.FunctionalTransitionSystem(solver)
    x = ts.make_state('x', bvsort4)
    mem = ts.make_state('mem', arrsort)

    u = c.Unroller(ts, solver)

    x0 = u.at_time(x, 0)
    assert x0 != x
    assert x0 == u.at_time(x, 0)

    assert x0 != u.at_time(x, 1)
    assert u.at_time(x, 1) == u.at_time(x, 1)

    constarr0 = solver.make_term(solver.make_term(0, bvsort8), arrsort)
    ts.constrain_init(solver.make_term(ss.primops.Equal,
                                       mem,
                                       constarr0))
    ts.assign_next(mem, solver.make_term(ss.primops.Store,
                                         mem,
                                         x,
                                         solver.make_term(1, bvsort8)))

    assert ts.init != u.at_time(ts.init, 0)
    assert u.at_time(ts.init, 0) == u.at_time(ts.init, 0)

    trans1 = u.at_time(ts.trans, 1)
    free_vars = get_free_vars(trans1)

    assert u.at_time(x, 1) in free_vars


@pytest.mark.parametrize("create_solver", ss.solvers.values())
def test_fts_unroller(create_solver):
    solver = create_solver()
    bvsort4 = solver.make_sort(ss.sortkinds.BV, 4)
    bvsort8 = solver.make_sort(ss.sortkinds.BV, 8)
    arrsort = solver.make_sort(ss.sortkinds.ARRAY, bvsort4, bvsort8)

    ts = c.RelationalTransitionSystem(solver)
    x = ts.make_state('x', bvsort4)
    mem = ts.make_state('mem', arrsort)

    u = c.Unroller(ts, solver)

    x0 = u.at_time(x, 0)
    assert x0 != x
    assert x0 == u.at_time(x, 0)

    assert x0 != u.at_time(x, 1)
    assert u.at_time(x, 1) == u.at_time(x, 1)

    constarr0 = solver.make_term(solver.make_term(0, bvsort8), arrsort)
    ts.constrain_init(solver.make_term(ss.primops.Equal,
                                       mem,
                                       constarr0))

    try:
        ts.assign_next(mem, solver.make_term(ss.primops.Store,
                                            mem,
                                            ts.next(x),
                                            solver.make_term(1, bvsort8)))
        assert False
    except:
        pass

    ts.constrain_trans(solver.make_term(ss.primops.Equal,
                                        ts.next(mem),
                                        solver.make_term(ss.primops.Store,
                                            mem,
                                            ts.next(x),
                                            solver.make_term(1, bvsort8))))

    assert ts.init != u.at_time(ts.init, 0)
    assert u.at_time(ts.init, 0) == u.at_time(ts.init, 0)

    trans1 = u.at_time(ts.trans, 1)
    free_vars = get_free_vars(trans1)

    assert u.at_time(x, 1) not in free_vars
    assert u.at_time(x, 2) in free_vars
