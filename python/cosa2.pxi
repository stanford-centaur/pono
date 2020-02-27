from cython.operator cimport dereference as dref, preincrement as inc
from libcpp.string cimport string
from libcpp.unordered_set cimport unordered_set
from libcpp.unordered_map cimport unordered_map
from libcpp.vector cimport vector

from cosa2 cimport TransitionSystem as c_TransitionSystem
from cosa2 cimport RelationalTransitionSystem as c_RelationalTransitionSystem
from cosa2 cimport FunctionalTransitionSystem as c_FunctionalTransitionSystem
from cosa2 cimport Property as c_Property
from cosa2 cimport Unroller as c_Unroller
from cosa2 cimport ProverResult as c_ProverResult
from cosa2 cimport UNKNOWN as c_UNKNOWN
from cosa2 cimport FALSE as c_FALSE
from cosa2 cimport TRUE as c_TRUE
from cosa2 cimport Prover as c_Prover
from cosa2 cimport Bmc as c_Bmc
from cosa2 cimport KInduction as c_KInduction
from cosa2 cimport BmcSimplePath as c_BmcSimplePath
from cosa2 cimport InterpolantMC as c_InterpolantMC

from smt_switch cimport SmtSolver, Sort, Term, c_Term, c_UnorderedTermMap

from enum import Enum

ctypedef const unordered_set[c_Term]* const_UnorderedTermSetPtr
ctypedef unordered_set[c_Term].const_iterator c_UnorderedTermSet_const_iterator

ctypedef const unordered_map[c_Term, c_Term]* const_UnorderedTermMapPtr
ctypedef unordered_map[c_Term, c_Term].const_iterator c_UnorderedTermMap_const_iterator

cdef class __AbstractTransitionSystem:
    cdef c_TransitionSystem* cts
    cdef SmtSolver _solver
    # Note: don't want to allow null TransitionSystems
    # means there's no way to instantiate a transition system without the solver
    def __cinit__(self, SmtSolver s):
        # if not specified, this creates a relational transition system under the hood
        self.cts = new c_RelationalTransitionSystem(s.css)
        self._solver = s

    def set_init(self, Term init):
        dref(self.cts).set_init(init.ct)

    def constrain_init(self, Term constraint):
        dref(self.cts).constrain_init(constraint.ct)

    def assign_next(self, Term state, Term val):
        dref(self.cts).assign_next(state.ct, val.ct)

    def add_invar(self, Term constraint):
        dref(self.cts).add_invar(constraint.ct)

    def constrain_inputs(self, Term constraint):
        dref(self.cts).constrain_inputs(constraint.ct)

    def add_constraint(self, Term constraint):
        dref(self.cts).add_constraint(constraint.ct)

    def name_term(self, str name, Term t):
        dref(self.cts).name_term(name.encode(), t.ct)

    def make_input(self, str name, Sort sort):
        cdef Term term = Term(self._solver)
        term.ct = dref(self.cts).make_input(name.encode(), sort.cs)
        return term

    def make_state(self, str name, Sort sort):
        cdef Term term = Term(self._solver)
        term.ct = dref(self.cts).make_state(name.encode(), sort.cs)
        return term

    def curr(self, Term t):
        cdef Term term = Term(self._solver)
        term.ct = dref(self.cts).curr(t.ct)
        return term

    def next(self, Term t):
        cdef Term term = Term(self._solver)
        term.ct = dref(self.cts).next(t.ct)
        return term

    def is_curr_var(self, Term sv):
        return dref(self.cts).is_curr_var(sv.ct)

    def is_next_var(self, Term sv):
        return dref(self.cts).is_next_var(sv.ct)

    @property
    def solver(self):
        return self._solver

    @property
    def states(self):
        states_set = set()

        cdef const_UnorderedTermSetPtr c_states_set = &dref(self.cts).states()
        cdef c_UnorderedTermSet_const_iterator it = c_states_set.const_begin()
        cdef c_UnorderedTermSet_const_iterator e  = c_states_set.const_end()

        cdef Term term
        while it != e:
            term = Term(self._solver)
            term.ct = dref(it)
            states_set.add(term)
            inc(it)

        return states_set

    @property
    def inputs(self):
        inputs_set = set()

        cdef const_UnorderedTermSetPtr c_inputs_set = &dref(self.cts).inputs()
        cdef c_UnorderedTermSet_const_iterator it = c_inputs_set.const_begin()
        cdef c_UnorderedTermSet_const_iterator e  = c_inputs_set.const_end()

        cdef Term term
        while it != e:
            term = Term(self._solver)
            term.ct = dref(it)
            inputs_set.add(term)
            inc(it)

        return inputs_set

    @property
    def state_updates(self):
        updates = {}

        cdef const_UnorderedTermMapPtr c_updates_map = &dref(self.cts).state_updates()
        cdef c_UnorderedTermMap_const_iterator it = c_updates_map.const_begin()
        cdef c_UnorderedTermMap_const_iterator e = c_updates_map.const_end()

        cdef Term k
        cdef Term v
        while it != e:
            k = Term(self._solver)
            v = Term(self._solver)
            k.ct = dref(it).first
            v.ct = dref(it).second
            updates[k] = v
            inc(it)

        return updates

    @property
    def named_terms(self):
        names2terms = {}

        cdef unordered_map[string, c_Term]* c_named_terms = &dref(self.cts).named_terms()
        cdef unordered_map[string, c_Term].const_iterator it = c_named_terms.const_begin()
        cdef unordered_map[string, c_Term].const_iterator e = c_named_terms.const_end()

        cdef Term term
        while it != e:
            term = Term(self._solver)
            term.ct = dref(it).second
            names2terms[(<string?> dref(it).first).decode()] = term
            inc(it)

        return names2terms

    @property
    def init(self):
        cdef Term term = Term(self._solver)
        term.ct = dref(self.cts).init()
        return term

    @property
    def trans(self):
        cdef Term term = Term(self._solver)
        term.ct = dref(self.cts).trans()
        return term

    def is_functional(self):
        return dref(self.cts).is_functional()


cdef class RelationalTransitionSystem(__AbstractTransitionSystem):
    def __cinit__(self, SmtSolver s):
        self.cts = new c_RelationalTransitionSystem(s.css)
        self._solver = s

    def set_behavior(self, Term init, Term trans):
        dref(<c_RelationalTransitionSystem * ?> self.cts).set_behavior(init.ct, trans.ct)

    def set_trans(self, Term trans):
        dref(<c_RelationalTransitionSystem * ?> self.cts).set_trans(trans.ct)

    def constrain_trans(self, Term constraint):
        dref(<c_RelationalTransitionSystem * ?> self.cts).constrain_trans(constraint.ct)


cdef class FunctionalTransitionSystem(__AbstractTransitionSystem):
    def __cinit__(self, SmtSolver s):
        self.cts = new c_FunctionalTransitionSystem(s.css)
        self._solver = s


cdef class Property:
    cdef c_Property* cp
    cdef __AbstractTransitionSystem ts
    def __cinit__(self, __AbstractTransitionSystem ts, Term p):
        self.cp = new c_Property(ts.cts[0], p.ct)
        self.ts = ts

    @property
    def prop(self):
        cdef Term p = Term(self.cts.solver)
        p.ct = dref(self.cp).prop()
        return p

    @property
    def transition_system(self):
        return self.ts


cdef class Unroller:
    cdef c_Unroller* cu
    cdef SmtSolver _solver
    def __cinit__(self, __AbstractTransitionSystem ts, SmtSolver s):
        self.cu = new c_Unroller(ts.cts[0], s.css)
        self._solver = s

    def at_time(self, Term t, unsigned int k):
        cdef Term term = Term(self._solver)
        term.ct = dref(self.cu).at_time(t.ct, k)
        return term

    def untime(self, Term t):
        cdef Term term = Term(self._solver)
        term.ct = dref(self.cu).untime(t.ct)
        return term


cdef class __AbstractProver:
    cdef c_Prover* cp
    cdef Property _property
    cdef SmtSolver _solver
    def __cinit__(self, Property p, SmtSolver s):
        self._property = p
        self._solver = s

    def initialize(self):
        dref(self.cp).initialize()

    def check_until(self, int k):
        '''
        Checks until bound k, returns True, False or None (if unknown)
        '''
        cdef int r = <int> dref(self.cp).check_until(k)
        if r == (<int> c_UNKNOWN):
            return None
        elif r == (<int> c_FALSE):
            return False
        elif r == (<int> c_TRUE):
            return True

    def witness(self):
        cdef vector[c_UnorderedTermMap] cw
        success = dref(self.cp).witness(cw)

        if not success:
            return None

        cdef Term kt
        cdef Term vt
        w = []
        for m in cw:
            d = dict()
            for elem in m:
                kt = Term(self._solver)
                kt.ct = (<c_Term?> elem.first)
                vt = Term(self._solver)
                vt.ct = (<c_Term?> elem.second)
                d[kt] = vt
            w.append(d)

        return w

    def prove(self):
        '''
        Tries to prove property unboundedly, returns True, False or None (if unknown)
        '''
        cdef int r = <int> dref(self.cp).prove()

        if r == (<int> c_UNKNOWN):
            return None
        elif r == (<int> c_FALSE):
            return False
        elif r == (<int> c_TRUE):
            return True

    @property
    def prop(self):
        return self._property


cdef class Bmc(__AbstractProver):
    def __cinit__(self, Property p, SmtSolver s):
        self.cp = new c_Bmc(p.cp[0], s.css)
        self._solver = s


cdef class KInduction(__AbstractProver):
    def __cinit__(self, Property p, SmtSolver s):
        self.cp = new c_KInduction(p.cp[0], s.css)
        self._solver = s


cdef class BmcSimplePath(__AbstractProver):
    def __cinit__(self, Property p, SmtSolver s):
        self.cp = new c_BmcSimplePath(p.cp[0], s.css)
        self._solver = s


cdef class InterpolantMC(__AbstractProver):
    def __cinit__(self, Property p, SmtSolver s, SmtSolver interp):
        self.cp = new c_InterpolantMC(p.cp[0], s.css, interp.css)
        self._solver = s
