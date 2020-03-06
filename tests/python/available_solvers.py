import smt_switch as ss

if hasattr(ss, 'create_msat_interpolator'):
    solver_and_interpolators = {'msat': (ss.create_msat_solver, ss.create_msat_interpolator)}
else:
    solver_and_interpolators = {}
bv_solvers = {name:fun for name, fun in ss.solvers.items()}
arith_solvers = {name:fun for name, fun in ss.solvers.items() if name != 'btor'}
