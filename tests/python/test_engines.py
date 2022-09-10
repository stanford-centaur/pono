import pytest
import smt_switch as ss
from smt_switch.sortkinds import BV
from smt_switch.primops import And, BVAdd, BVSub, Equal, Ite
import pono
import available_solvers

from typing import Tuple

def build_simple_alu_fts(s:ss.SmtSolver)->Tuple[pono.Property, pono.FunctionalTransitionSystem]:
    '''
    Creates a simple alu transition system
    @param s - an SmtSolver from smt_switch
    @return a property
    '''

    # Instantiate a functional transition system
    fts = pono.FunctionalTransitionSystem(s)

    # Create a bit-vector sorts
    bvsort1 = s.make_sort(BV, 1)
    bvsort8 = s.make_sort(BV, 8)

    # Create the states
    cfg = fts.make_statevar('cfg', bvsort1)
    spec_res = fts.make_statevar('spec_res', bvsort8)
    imp_res  = fts.make_statevar('imp_res', bvsort8)

    # Create the inputs
    a = fts.make_inputvar('a', bvsort8)
    b = fts.make_inputvar('b', bvsort8)

    # Add logic for cfg
    ## Start at 0
    fts.constrain_init(s.make_term(Equal, cfg, s.make_term(0, bvsort1)))
    ## Keeps the same value
    fts.assign_next(cfg, cfg)

    # Set logic for results
    ## they start equal
    fts.constrain_init(s.make_term(Equal, spec_res, imp_res))
    ## spec_res is the sum: spec_res' = a + b
    fts.assign_next(spec_res, s.make_term(BVAdd, a, b))
    ## depends on the configuration: imp_res' == (cfg == 0) ? a + b : a - b
    fts.assign_next(imp_res, s.make_term(Ite,
                                     s.make_term(Equal, cfg, s.make_term(0, bvsort1)),
                                     s.make_term(BVAdd, a, b),
                                     s.make_term(BVSub, a, b)))

    # Create a property: (spec_cnt == imp_cnt - 1)
    prop = pono.Property(s, s.make_term(Equal, spec_res, imp_res))
    return prop, fts

@pytest.mark.parametrize("create_solver", ss.solvers.values())
def test_bmc(create_solver):
    s = create_solver(False)
    s.set_opt('produce-models', 'true')
    s.set_opt('incremental', 'true')
    prop, ts = build_simple_alu_fts(s)

    bmc = pono.Bmc(prop, ts, s)
    res = bmc.check_until(10)

    assert res is None, "BMC shouldn't be able to solve"

@pytest.mark.parametrize("create_solver", ss.solvers.values())
def test_kind(create_solver):
    s = create_solver(False)
    s.set_opt('produce-models', 'true')
    s.set_opt('incremental', 'true')
    prop, ts = build_simple_alu_fts(s)

    kind = pono.KInduction(prop, ts, s)
    res = kind.check_until(1)

    # NOTE: k-induction uses a stronger termination criterion based on
    # initial states by default (newly added), which allows it to prove the
    # property with only one unrolling.
    assert res is True, "KInduction should be able to solve this property"


@pytest.mark.parametrize("solver_and_interpolator", available_solvers.solver_and_interpolators.values())
def test_interp(solver_and_interpolator):
    s = solver_and_interpolator[0](False)
    s.set_opt('produce-models', 'true')
    s.set_opt('incremental', 'true')

    prop, ts = build_simple_alu_fts(s)

    ic3ia = pono.IC3IA(prop, ts, s)
    res = ic3ia.check_until(10)

    assert res is True, "IC3IA be able to solve this property"


@pytest.mark.parametrize("solver_and_interpolator", available_solvers.solver_and_interpolators.values())
def test_interp(solver_and_interpolator):
    s = solver_and_interpolator[0](False)
    s.set_opt('produce-models', 'true')
    s.set_opt('incremental', 'true')

    prop, ts = build_simple_alu_fts(s)

    interp = pono.InterpolantMC(prop, ts, s)
    res = interp.check_until(10)

    assert res is True, "InterpolantMC be able to solve this property"

@pytest.mark.parametrize("create_solver", ss.solvers.values())
def test_kind_inductive_prop(create_solver):
    s = create_solver(False)
    s.set_opt('produce-models', 'true')
    s.set_opt('incremental', 'true')
    prop, ts = build_simple_alu_fts(s)

    states = {str(sv):sv for sv in ts.statevars}

    prop = pono.Property(s, s.make_term(And,
                                        s.make_term(Equal, states['cfg'],
                                                    s.make_term(0, s.make_sort(BV, 1))),
                                        prop.prop))

    kind = pono.KInduction(prop, ts, s)
    res = kind.check_until(10)

    assert res is True, "KInduction should be able to solve this manually strengthened property"
