from cython.operator cimport dereference as dref, preincrement as inc
from libcpp.string cimport string
from libcpp.unordered_set cimport unordered_set
from libcpp.unordered_map cimport unordered_map

from cosa2 cimport TransitionSystem as c_TransitionSystem
from cosa2 cimport FunctionalTransitionSystem as c_FunctionalTransitionSystem
from cosa2 cimport Property as c_Property

from smt_switch cimport SmtSolver, Sort, Term, c_Term, c_UnorderedTermMap

ctypedef const unordered_set[c_Term]* const_UnorderedTermSetPtr
ctypedef unordered_set[c_Term].const_iterator c_UnorderedTermSet_const_iterator

ctypedef const unordered_map[c_Term, c_Term]* const_UnorderedTermMapPtr
ctypedef unordered_map[c_Term, c_Term].const_iterator c_UnorderedTermMap_const_iterator

cdef class TransitionSystem:
    cdef c_TransitionSystem* cts
    cdef SmtSolver _solver
    # Note: don't want to allow null TransitionSystems
    # means there's no way to instantiate a transition system without the solver
    def __cinit__(self, SmtSolver s):
        self.cts = new c_TransitionSystem(s.css)
        self._solver = s

    def set_behavior(self, Term init, Term trans):
        dref(self.cts).set_behavior(init.ct, trans.ct)

    def set_init(self, Term init):
        dref(self.cts).set_init(init.ct)

    def constrain_init(self, Term constraint):
        dref(self.cts).constrain_init(constraint.ct)

    def set_trans(self, Term trans):
        dref(self.cts).set_trans(trans.ct)

    def constrain_trans(self, Term constraint):
        dref(self.cts).constrain_trans(constraint.ct)

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


cdef class FunctionalTransitionSystem(TransitionSystem):
    def __cinit__(self, SmtSolver s):
        self.cts = new c_FunctionalTransitionSystem(s.css)
        self._solver = s


cdef class Property:
    cdef c_Property* cp
    cdef TransitionSystem ts
    def __cinit__(self, TransitionSystem ts, Term p):
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
