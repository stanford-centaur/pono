import pytest
import smt_switch as ss
import pono


@pytest.mark.parametrize("create_solver", ss.solvers.values())
def test_cons_fts(create_solver):
    solver = create_solver(False)
    boolsort = solver.make_sort(ss.sortkinds.BOOL)
    bvsort8 = solver.make_sort(ss.sortkinds.BV, 8)

    fts = pono.FunctionalTransitionSystem(solver);

    a = fts.make_inputvar("a", bvsort8);
    b = fts.make_inputvar("b", bvsort8);
    c = fts.make_inputvar("c", bvsort8);
    d = fts.make_inputvar("d", boolsort);

    regres = fts.make_statevar("regres", bvsort8);
    counter = fts.make_statevar("counter", bvsort8);

    zero = fts.make_term(0, bvsort8);
    one = fts.make_term(1, bvsort8);
    fts.constrain_init(fts.make_term(ss.primops.Equal, counter, zero));
    fts.assign_next(counter, fts.make_term(ss.primops.BVAdd, counter, one));
    fts.assign_next(regres, fts.make_term(ss.primops.BVAdd, a, b));
    fts.add_constraint(fts.make_term(ss.primops.Equal, d, fts.make_term(ss.primops.BVUle, a, b)));

    out = fts.make_term(ss.primops.BVSub, regres, one);
    fts.name_term("out", out);

    pono.coi_reduction(fts, [ out ])

    assert regres in fts.statevars
    assert len(fts.inputvars) == 3
    assert a in fts.inputvars
    assert b in fts.inputvars
    assert c not in fts.inputvars
    assert d in fts.inputvars
    assert "a" in fts.named_terms
    assert "c" not in fts.named_terms
