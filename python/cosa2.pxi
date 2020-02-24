from cosa2 cimport TransitionSystem as c_TransitionSystem

from smt_switch import SmtSolver, Sort, Term


cdef class TransitionSystem:
    cdef c_TransitionSystem cts
    def __cinit__(self, SmtSolver s):
        cts = c_TransitionSystem(s.css)

    def set_behavior(self, Term init, Term trans):
        self.cts.set_behavior(init.ct, trans.ct)

    def set_init(self, Term init):
        self.cts.set_init(init.ct)

    def constrain_init(self, Term constraint):
        self.cts.constrain_init(constraint.ct)

    def set_trans(self, Term trans):
        self.cts.set_trans(trans.ct)

    def constrain_trans(self, Term constraint):
        self.cts.constrain_trans(constraint.ct)

    def assign_next(self, Term state, Term val):
        self.cts.assign_next(state.ct, val.ct)

    def add_invar(self, Term constraint):
        self.cts.add_invar(constraint.ct)

    def constrain_inputs(self, Term constraint):
        self.cts.constrain_inputs(constraint.ct)

    def add_constraint(self, Term constraint):
        self.cts.add_constraint(constraint.ct)

    def name_term(self, str name, Term t):
        self.cts.name_term(name.encode(), t.ct)

    def make_input(self, str name, Sort sort):
        cdef Term term = Term()
        term.ct = self.cts.make_input(name.encode(), sort.cs)
        return term

    def make_state(self, str name, Sort sort):
        cdef Term term = Term()
        term.ct = self.cts.make_state(name.encode(), sort.cs)
        return term

    def curr(self, Term t):
        cdef Term term = Term()
        term.ct = self.cts.curr(t.ct)
        return term

    def next(self, Term t):
        cdef Term term = Term()
        term.ct = self.cts.next(t.ct)
        return term

    def is_curr_var(self, Term sv):
        return self.cts.is_curr_var(sv.ct)

    def is_next_var(self, Term sv):
        return self.cts.is_next_var(sv.ct)

    @property
    def solver(self):
        cdef SmtSolver ss = SmtSolver()
        ss.css = self.cts.solver()
        return ss

    # TODO: uncomment these (might need more iteration operators)
    # @property
    # def states(self):
    #     states_set = set()

    #     cdef Term term = Term()
    #     for s in self.cts.states():
    #         term.ct = s
    #         states_set.insert(term)

    #     return states_set

    # @property
    # def inputs(self):
    #     inputs_set = set()

    #     cdef Term term = Term()
    #     for s in self.cts.inputs():
    #         term.ct = s
    #         inputs_set.insert(term)

    #     return inputs_set

    # @property
    # def state_updates(self):
    #     updates = {}
    #     cdef Term k = Term()
    #     cdef Term v = Term()

    #     for elem in self.cts.state_updates():
    #         k.ct = elem.first
    #         v.ct = elem.second
    #         updates[k] = v

    #     return updates

    # @property
    # def named_terms(self):
    #     names2terms = {}

    #     cdef Term term = Term()
    #     for elem in self.cts.named_terms():
    #         term.ct = elem.second
    #         names2terms[elem.first.decode()] = term

    #     return names2terms

    @property
    def init(self):
        cdef Term term = Term()
        term.ct = self.cts.init()
        return term

    @property
    def trans(self):
        cdef Term term = Term()
        term.ct = self.cts.trans()
        return term

    def is_functional(self):
        return self.cts.is_functional()
