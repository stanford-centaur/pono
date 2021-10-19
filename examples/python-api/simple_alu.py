#!/usr/bin/env python3
import argparse
import pono
import smt_switch as ss
from smt_switch.primops import And, BVAdd, BVSub, Equal, Ite
from smt_switch.sortkinds import BOOL, BV

def build_simple_alu_fts(s:ss.SmtSolver)->pono.Property:
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
    ## imp_res depends on the configuration: imp_res' == (cfg == 0) ? a + b : a - b
    fts.assign_next(imp_res, s.make_term(Ite,
                                     s.make_term(Equal, cfg, s.make_term(0, bvsort1)),
                                     s.make_term(BVAdd, a, b),
                                     s.make_term(BVSub, a, b)))

    # Create a property: spec_res == imp_res
    prop = pono.Property(s, s.make_term(Equal,
                                        spec_res,
                                        imp_res))
    return prop, fts


def k_induction_attempt():
    # Create an smt_switch.SmtSolver with Boolector as the backend
    # and no logging
    s = ss.create_btor_solver(False)
    s.set_opt('produce-models', 'true')
    s.set_opt('incremental', 'true')
    prop, fts = build_simple_alu_fts(s)

    print('\n============== Running k-induction ==============')
    print('INIT\n\t{}'.format(fts.init))
    print('TRANS\n\t{}'.format(fts.trans))
    print('PROP\n\t{}'.format(prop.prop))

    # Create KInduction engine -- using same solver (in future can change the solver)
    kind = pono.KInduction(prop, fts, s)
    res = kind.check_until(20)

    print(res)

    assert res is None, "Expecting k-induction not to prove property in 20 steps"
    print("KInduction returned unknown")


def interpolant_attempt():
    # Create solver and interpolator using MathSAT
    # and no logging for the solver
    s = ss.create_msat_solver(False)
    s.set_opt('produce-models', 'true')
    s.set_opt('incremental', 'true')
    prop, fts = build_simple_alu_fts(s)

    print('\n============== Running Interpolant-based Model Checking ==============')
    print('INIT\n\t{}'.format(fts.init))
    print('TRANS\n\t{}'.format(fts.trans))
    print('PROP\n\t{}'.format(prop.prop))

    # Create InterpolantMC engine
    itpmc = pono.InterpolantMC(prop, fts, s)
    res = itpmc.check_until(20)

    print(res)

    assert res is True, "Expecting InterpolantMC to prove the property"
    print("InterpolantMC returned true")

def k_induction_attempt_inductive():
    # Create an smt_switch.SmtSolver with Boolector as the backend
    # and no logging
    s = ss.create_btor_solver(False)
    s.set_opt('produce-models', 'true')
    s.set_opt('incremental', 'true')
    prop, fts = build_simple_alu_fts(s)

    # store sets of states in a dictionary for accessing below
    states = {str(sv):sv for sv in fts.statevars}

    # make the property inductive manually
    prop = pono.Property(s,
                          s.make_term(And,
                                      s.make_term(Equal,
                                                  states['cfg'],
                                                  s.make_term(0, s.make_sort(BV, 1))),
                                      prop.prop))

    print('\n============== Running k-induction on inductively strengthened property ==============')
    print('INIT\n\t{}'.format(fts.init))
    print('TRANS\n\t{}'.format(fts.trans))
    print('PROP\n\t{}'.format(prop.prop))

    # Create KInduction engine -- using same solver (in future can change the solver)
    kind = pono.KInduction(prop, fts, s)
    res = kind.check_until(20)

    print(res)

    assert res is True, "Expecting k-induction to prove the inductively strengthened property"
    print("KInduction returned true")


approaches = {
    'kind': k_induction_attempt,
    'interp': interpolant_attempt,
    'kind-manual': k_induction_attempt_inductive
}

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Pono SimpleALU example')
    parser.add_argument('approach', choices=['kind', 'interp', 'kind-manual'],
                        help='Select the approach: k-induction, interpolant-based,'
                        ' or k-induction with a manually strengthened property')
    parser.add_argument('-v', '--verbosity', type=int, default=0)
    args = parser.parse_args()

    pono.set_global_logger_verbosity(args.verbosity)

    approaches[args.approach]()
