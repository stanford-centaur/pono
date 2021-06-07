###############################################################
# \file test_pseudo_init_and_prop.py
# \verbatim
# Top contributors (to current version):
#   Makai Mann
# This file is part of the smt-switch project.
# Copyright (c) 2021 by the authors listed in the file AUTHORS
# in the top-level source directory) and their institutional affiliations.
# All rights reserved.  See the file LICENSE in the top-level source
# directory for licensing information.\endverbatim
#
# \brief Test the pseudo_init_and_prop modification function.
#        This method is general but is most useful for improving
#        IC3IA performance. This test is based on the C++ tests.
#

import pytest
import smt_switch as ss
import pono
from typing import Set


@pytest.mark.parametrize("create_solver", ss.solvers.values())
def test_counter_system_safe(create_solver):
    solver = create_solver(False)
    solver.set_opt("produce-models", "true")
    solver.set_opt("incremental", "true")
    fts = pono.FunctionalTransitionSystem(solver)
    bvsort8 = fts.make_sort(ss.sortkinds.BV, 8)
    x = fts.make_statevar("x", bvsort8)
    one = fts.make_term(1, bvsort8)
    ten = fts.make_term(10, bvsort8)
    fts.assign_next(x, solver.make_term(ss.primops.Ite,
                                        solver.make_term(ss.primops.BVUlt,
                                                         x,
                                                         ten),
                                        solver.make_term(ss.primops.BVAdd,
                                                         x,
                                                         one),
                                        x))
    fts.set_init(solver.make_term(ss.primops.Equal,
                                  x,
                                  solver.make_term(0, bvsort8)))
    prop_term = fts.make_term(ss.primops.BVUle, x, ten)
    rts_new = pono.pseudo_init_and_prop(fts, prop_term)
    p = pono.Property(solver, prop_term)

    kind = pono.KInduction(p, rts_new, solver)
    res = kind.check_until(11)
    assert res


@pytest.mark.parametrize("create_solver", ss.solvers.values())
def test_trivial_unsafe(create_solver):
    solver = create_solver(False)
    solver.set_opt("produce-models", "true")
    solver.set_opt("incremental", "true")
    rts = pono.RelationalTransitionSystem(solver)
    bvsort8 = rts.make_sort(ss.sortkinds.BV, 8)
    x = rts.make_statevar('x', bvsort8)
    prop_term = solver.make_term(ss.primops.BVUle, x, solver.make_term(10, bvsort8))

    rts.set_trans(rts.make_term(False))


    rts_new = pono.pseudo_init_and_prop(rts, prop_term)
    p = pono.Property(solver, prop_term)

    kind = pono.KInduction(p, rts_new, solver)
    res = kind.check_until(11)
    assert not res
