from cython.operator cimport dereference as dref, preincrement as inc
from libc.stdint cimport uintptr_t
from libcpp.string cimport string
from libcpp.unordered_set cimport unordered_set
from libcpp.unordered_map cimport unordered_map
from libcpp.vector cimport vector

from pono_imp cimport TransitionSystem as c_TransitionSystem
from pono_imp cimport RelationalTransitionSystem as c_RelationalTransitionSystem
from pono_imp cimport FunctionalTransitionSystem as c_FunctionalTransitionSystem
from pono_imp cimport Property as c_Property
from pono_imp cimport Unroller as c_Unroller
from pono_imp cimport ProverResult as c_ProverResult
from pono_imp cimport UNKNOWN as c_UNKNOWN
from pono_imp cimport FALSE as c_FALSE
from pono_imp cimport TRUE as c_TRUE
from pono_imp cimport Prover as c_Prover
from pono_imp cimport Bmc as c_Bmc
from pono_imp cimport KInduction as c_KInduction
from pono_imp cimport BmcSimplePath as c_BmcSimplePath
from pono_imp cimport InterpolantMC as c_InterpolantMC
from pono_imp cimport BTOR2Encoder as c_BTOR2Encoder
IF WITH_COREIR == "ON":
    from pono_imp cimport Module as c_Module
    from pono_imp cimport CoreIREncoder as c_CoreIREncoder
from pono_imp cimport set_global_logger_verbosity as c_set_global_logger_verbosity
from pono_imp cimport get_free_symbols as c_get_free_symbols

from smt_switch cimport SmtSolver, Sort, Term, c_Term, c_UnorderedTermMap

from enum import Enum

PYCOREIR_AVAILABLE=False
try:
    import coreir
    import ctypes
    PYCOREIR_AVAILABLE=True
except:
    print("Warning: Pono built with CoreIR support but coreir python module not found")

ctypedef unordered_set[c_Term] c_UnorderedTermSet
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

    def make_inputvar(self, str name, Sort sort):
        cdef Term term = Term(self._solver)
        term.ct = dref(self.cts).make_inputvar(name.encode(), sort.cs)
        return term

    def make_statevar(self, str name, Sort sort):
        cdef Term term = Term(self._solver)
        term.ct = dref(self.cts).make_statevar(name.encode(), sort.cs)
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
    def statevars(self):
        states_set = set()

        cdef const_UnorderedTermSetPtr c_states_set = &dref(self.cts).statevars()
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
    def inputvars(self):
        inputs_set = set()

        cdef const_UnorderedTermSetPtr c_inputs_set = &dref(self.cts).inputvars()
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
        cdef Term p = Term(self.ts.solver)
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


cdef class BTOR2Encoder:
    cdef c_BTOR2Encoder * cbe
    def __cinit__(self, str filename, __AbstractTransitionSystem ts):
        self.cbe = new c_BTOR2Encoder(filename.encode(), dref(ts.cts))

IF WITH_COREIR == "ON":
    cdef class CoreIREncoder:
        cdef c_CoreIREncoder * cbe
        def __cinit__(self, mod, RelationalTransitionSystem ts):
            cdef uintptr_t adr
            if isinstance(mod, str):
                self.cbe = new c_CoreIREncoder((<string?> (mod.encode())), dref((<c_RelationalTransitionSystem *> ts.cts)))
            elif hasattr(mod, "ptr"):
                adr = <uintptr_t> ctypes.addressof(mod.ptr.contents)
                self.cbe = new c_CoreIREncoder((<c_Module *> adr), dref((<c_RelationalTransitionSystem *> ts.cts)))
            else:
                raise ValueError("CoreIR encoder takes a pycoreir Context or a filename but got {}".format(mod))



def set_global_logger_verbosity(int v):
    c_set_global_logger_verbosity(v)


def get_free_symbols(Term term):
    cdef c_UnorderedTermSet out_symbols
    c_get_free_symbols(term.ct, out_symbols)

    python_out_set = set()
    for s in out_symbols:
        t = Term(term._solver)
        t.ct = s
        python_out_set.add(t)

    return python_out_set
