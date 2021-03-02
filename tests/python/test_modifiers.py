import pytest
import smt_switch as ss
import pono

@pytest.mark.parametrize("create_solver", ss.solvers.values())
def test_history_modifier(create_solver):
    solver = create_solver(False)
    solver.set_opt("incremental", "true")
    solver.set_opt("produce-models", "true")
    bvsort8 = solver.make_sort(ss.sortkinds.BV, 8)

    ts = pono.FunctionalTransitionSystem(solver)
    counter = ts.make_statevar("counter", bvsort8)
    ts.constrain_init(ts.make_term(ss.primops.Equal, counter, ts.make_term(0, bvsort8)))
    ts.assign_next(counter, ts.make_term(ss.primops.BVAdd, counter, ts.make_term(1, bvsort8)))

    hm = pono.HistoryModifier(ts)
    counter_delay_2 = hm.get_hist(counter, 2)

    p = pono.Property(solver, ts.make_term(ss.primops.BVUlt, counter, ts.make_term(5, bvsort8)))

    bmc = pono.Bmc(p, ts, solver)
    res = bmc.check_until(10)
    assert not res, "should be false"

    witness = bmc.witness()

    checked_at_least_one = False
    for i, m in enumerate(witness):
        if i > 1:
            checked_at_least_one = True
            # checking semantics of history modifier delay
            assert witness[i][counter_delay_2] == witness[i-2][counter]

    assert checked_at_least_one
