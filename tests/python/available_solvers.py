# Copyright (c) 2020 by the authors listed in the file AUTHORS in the top-level
# source directory and their institutional affiliations. See the file LICENSE
# in the top-level source directory for licensing information.
#
# This file is part of the pono project.
import smt_switch as ss

if hasattr(ss, "create_msat_interpolator"):
    solver_and_interpolators = {
        "msat": (ss.create_msat_solver, ss.create_msat_interpolator)
    }
else:
    solver_and_interpolators = {}
bv_solvers = dict(ss.solvers.items())
arith_solvers = {name: fun for name, fun in ss.solvers.items() if name != "btor"}
