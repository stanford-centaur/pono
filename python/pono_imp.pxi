from cython.operator cimport dereference as dref, preincrement as inc
from libc.stdint cimport uintptr_t
from libcpp cimport bool as cbool
from libcpp.pair cimport pair
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
from pono_imp cimport IC3 as c_IC3
from pono_imp cimport IC3Bits as c_IC3Bits
from pono_imp cimport IC3IA as c_IC3IA
from pono_imp cimport IC3SA as c_IC3SA
from pono_imp cimport InterpolantMC as c_InterpolantMC
from pono_imp cimport ModelBasedIC3 as c_ModelBasedIC3
IF WITH_MSAT_IC3IA == "ON":
    from pono_imp cimport MsatIC3IA as c_MsatIC3IA
from pono_imp cimport BTOR2Encoder as c_BTOR2Encoder
from pono_imp cimport SMVEncoder as c_SMVEncoder
from pono_imp cimport VMTEncoder as c_VMTEncoder
IF WITH_COREIR == "ON":
    from pono_imp cimport Module as c_Module
    from pono_imp cimport CoreIREncoder as c_CoreIREncoder
from pono_imp cimport PonoOptions as c_PonoOptions
from pono_imp cimport HistoryModifier as c_HistoryModifier
from pono_imp cimport StaticConeOfInfluence as c_StaticConeOfInfluence
from pono_imp cimport add_prop_monitor as c_add_prop_monitor
from pono_imp cimport VCDWitnessPrinter as c_VCDWitnessPrinter
from pono_imp cimport pseudo_init_and_prop as c_pseudo_init_and_prop
from pono_imp cimport prop_in_trans as c_prop_in_trans
from pono_imp cimport set_global_logger_verbosity as c_set_global_logger_verbosity
from pono_imp cimport check_invar as c_check_invar

from smt_switch cimport SmtSolver, PrimOp, Op, c_SortKind, SortKind, \
    c_Sort, c_SortVec, Sort, Term, c_Term, c_TermVec, c_UnorderedTermMap

from enum import Enum

PYCOREIR_AVAILABLE=False
IF WITH_COREIR == "ON":
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
    # this pointer is allocated and deallocated by derived classes
    cdef c_TransitionSystem* cts
    cdef SmtSolver _solver
    def __cinit__(self, SmtSolver s):
        pass

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

    def add_constraint(self, Term constraint, bint to_init_and_next=True):
        dref(self.cts).add_constraint(constraint.ct, to_init_and_next)

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

    def is_input_var(self, Term iv):
        return dref(self.cts).is_input_var(iv.ct)

    def get_name(self, Term t):
        return dref(self.cts).get_name(t.ct).decode()

    def lookup(self, str name):
        cdef Term term = Term(self._solver)
        term.ct = dref(self.cts).lookup(name.encode())
        return term

    def add_statevar(self, Term cv, Term nv):
        dref(self.cts).add_statevar(cv.ct, nv.ct)

    def add_inputvar(self, Term v):
        dref(self.cts).add_inputvar(v.ct)

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
    def constraints(self):
        convec = []
        cdef vector[pair[c_Term, cbool]] c_cons = dref(self.cts).constraints()
        for e in c_cons:
            python_term = Term(self._solver)
            python_term.ct = e.first
            convec.append((python_term, e.second))
        return convec

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

    def is_deterministic(self):
        return dref(self.cts).is_deterministic()

    def drop_state_updates(self, list svs):
        '''
        EXPERTS ONLY
        Drop the state updates for these state variables
        '''
        cdef c_TermVec c_svs
        c_svs.reserve(len(svs))
        for sv in svs:
            c_svs.push_back((<Term?> sv).ct)

        dref(self.cts).drop_state_updates(c_svs)

    def promote_inputvar(self, Term iv):
        '''
        EXPERTS ONLY
        Turns an input variable into a state variable
          IMPORTANT: this does not retroactively change constraints
          e.g. if a constraint was not added to init because it
          contains an input variable

        @param iv the input variable to promote

        The input variable iv stays the same, but it will now
          be registered as a state variable.
        '''
        dref(self.cts).promote_inputvar(iv.ct)

    def replace_terms(self, dict to_replace):
        '''
        EXPERTS ONLY
        Use the provided dictionary to substitute all the terms with the mapping in the transition system
        Can be used to cut out pieces of the transition system that are unneeded. Note, there are no guarantees,
        the user is responsible for maintaining any semantics of the system that they want.

        Modifies the transition system in place.
        '''

        cdef c_UnorderedTermMap utm

        for k, v in to_replace.items():
            utm[(<Term?> k).ct] = (<Term?> v).ct

        dref(self.cts).replace_terms(utm)

    def make_sort(self, arg0, arg1=None, arg2=None, arg3=None):
        cdef Sort s = Sort(self._solver)
        cdef c_SortKind sk
        cdef c_SortVec csv

        if isinstance(arg0, str):
            s.cs = dref(self.cts).make_sort(<const string> arg0.encode(), <int?> arg1)
        elif isinstance(arg0, SortKind):
            sk = (<SortKind> arg0).sk
            if arg1 is None:
                s.cs = dref(self.cts).make_sort(sk)
            elif isinstance(arg1, int) and arg2 is None:
                s.cs = dref(self.cts).make_sort(sk, <int> arg1)
            elif isinstance(arg1, Sort) and arg2 is None:
                s.cs = dref(self.cts).make_sort(sk, (<Sort> arg1).cs)
            elif isinstance(arg1, list) and arg2 is None:
                for a in arg1:
                    csv.push_back((<Sort?> a).cs)
                s.cs = dref(self.cts).make_sort(sk, csv)
            elif arg3 is None:
                s.cs = dref(self.cts).make_sort(sk, (<Sort?> arg1).cs, (<Sort?> arg2).cs)
            elif arg3 is not None:
                s.cs = dref(self.cts).make_sort(sk, (<Sort?> arg1).cs,
                                                    (<Sort?> arg2).cs,
                                                    (<Sort?> arg3).cs)
            else:
                raise ValueError("Cannot find matching function for {}".format([type(a)
                                                                                for a in
                                                                                [arg0, arg1, arg2, arg3]]))
        else:
            raise ValueError("Cannot find matching function for {}".format([type(a)
                                                                            for a in
                                                                            [arg0, arg1, arg2, arg3]]))
        return s

    def make_term(self, op_or_val, *args):
        cdef Term term = Term(self._solver)
        cdef c_TermVec ctv

        if isinstance(op_or_val, PrimOp):
            op_or_val = Op(op_or_val)

        # expand a list argument
        if len(args) > 0:
            if (isinstance(args[0], list) and len(args) > 1) or \
               any([isinstance(a, list) for a in args[1:]]):
                raise ValueError("Cannot call make_term with signature {}".format([type(a) for a in args]))
            elif isinstance(args[0], list):
                # expand arguments in list to be args
                args = args[0]

        if isinstance(op_or_val, Op):
            if not op_or_val:
                raise ValueError("Got a null Op in make_term")

            if not args:
                raise ValueError("Can't call make_term with an Op ({}) and no arguments".format(op_or_val))

            for a in args:
                ctv.push_back((<Term?> a).ct)
            term.ct = dref(self.cts).make_term((<Op> op_or_val).op, ctv)
        elif isinstance(op_or_val, bool) and len(args) == 0:
            term.ct = dref(self.cts).make_term(<bint> op_or_val)
        elif isinstance(op_or_val, str) and len(args) == 1 and isinstance(args[0], Sort):
            term.ct = dref(self.cts).make_term(<const string> op_or_val.encode(), (<Sort> args[0]).cs)
        elif isinstance(op_or_val, str) and len(args) == 2 and isinstance(args[0], Sort):
            term.ct = dref(self.cts).make_term(<const string> op_or_val.encode(),
                                               (<Sort> args[0]).cs,
                                               <int?> args[1])
        elif isinstance(op_or_val, int) and len(args) == 1 and isinstance(args[0], Sort):
            # always use the string representation of integers (to handle large ints)
            term.ct = dref(self.cts).make_term((<const string?> str(op_or_val).encode()), (<Sort> args[0]).cs)
        elif isinstance(op_or_val, Term) and len(args) == 1 and isinstance(args[0], Sort):
            # this is for creating a constant array
            term.ct = dref(self.cts).make_term((<Term?> op_or_val).ct, (<Sort?> args[0]).cs)
        else:
            raise ValueError("Couldn't find matching function for {}".format([type(a)
                                                                              for a in [op_or_val] + args]))
        return term


cdef class RelationalTransitionSystem(__AbstractTransitionSystem):
    def __cinit__(self, SmtSolver s):
        self.cts = new c_RelationalTransitionSystem(s.css)
        self._solver = s

    def __dealloc__(self):
        del self.cts

    def set_behavior(self, Term init, Term trans):
        dref(<c_RelationalTransitionSystem * ?> self.cts).set_behavior(init.ct, trans.ct)

    def set_trans(self, Term trans):
        dref(<c_RelationalTransitionSystem * ?> self.cts).set_trans(trans.ct)

    def constrain_trans(self, Term constraint):
        dref(<c_RelationalTransitionSystem * ?> self.cts).constrain_trans(constraint.ct)

    def __deepcopy__(self, memo):
        cdef c_RelationalTransitionSystem * c_res = new c_RelationalTransitionSystem(dref(self.cts))
        res = RelationalTransitionSystem(self._solver)
        del res.cts
        res.cts = c_res
        return res


cdef class FunctionalTransitionSystem(__AbstractTransitionSystem):
    def __cinit__(self, SmtSolver s):
        self.cts = new c_FunctionalTransitionSystem(s.css)
        self._solver = s

    def __dealloc__(self):
        del self.cts

    def __deepcopy__(self, memo):
        cdef c_FunctionalTransitionSystem * c_res = new c_FunctionalTransitionSystem(dref(self.cts))
        res = FunctionalTransitionSystem(self._solver)
        del res.cts
        res.cts = c_res
        return res


cdef class Property:
    cdef c_Property* cp
    cdef SmtSolver _solver
    def __cinit__(self, SmtSolver s, Term p):
        self.cp = new c_Property(s.css, p.ct)
        self._solver = s

    def __dealloc__(self):
        del self.cp

    @property
    def prop(self):
        cdef Term p = Term(self._solver)
        p.ct = dref(self.cp).prop()
        return p

    @property
    def solver(self):
        return self._solver


cdef class Unroller:
    cdef c_Unroller* cu
    cdef SmtSolver _solver
    def __cinit__(self, __AbstractTransitionSystem ts, str time_id='@'):
        self.cu = new c_Unroller(ts.cts[0], time_id.encode())
        self._solver = ts._solver

    def __dealloc__(self):
        del self.cu

    def at_time(self, Term t, unsigned int k):
        cdef Term term = Term(self._solver)
        term.ct = dref(self.cu).at_time(t.ct, k)
        return term

    def untime(self, Term t):
        cdef Term term = Term(self._solver)
        term.ct = dref(self.cu).untime(t.ct)
        return term

    def get_var_time(self, Term v):
        return dref(self.cu).get_var_time(v.ct)


cdef class __AbstractProver:
    # this pointer is allocated and deallocated by derived classes
    cdef c_Prover* cp
    cdef Property _property
    cdef __AbstractTransitionSystem _ts
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

    def invar(self):
        cdef Term inv = Term(self._solver)
        inv.ct = dref(self.cp).invar()
        return inv

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
    def __cinit__(self, Property p, __AbstractTransitionSystem ts, SmtSolver s):
        self.cp = new c_Bmc(p.cp[0], ts.cts[0], s.css)
        self._ts = ts
        self._solver = s

    def __dealloc__(self):
        del self.cp


cdef class KInduction(__AbstractProver):
    def __cinit__(self, Property p, __AbstractTransitionSystem ts, SmtSolver s):
        self.cp = new c_KInduction(p.cp[0], ts.cts[0], s.css)
        self._ts = ts
        self._solver = s

    def __dealloc__(self):
        del self.cp


cdef class BmcSimplePath(__AbstractProver):
    def __cinit__(self, Property p, __AbstractTransitionSystem ts, SmtSolver s):
        self.cp = new c_BmcSimplePath(p.cp[0], ts.cts[0], s.css)
        self._ts = ts
        self._solver = s

    def __dealloc__(self):
        del self.cp


cdef class IC3(__AbstractProver):
    '''
    Boolean IC3 variant
    '''
    def __cinit__(self, Property p, __AbstractTransitionSystem ts, SmtSolver s):
        self.cp = new c_IC3(p.cp[0], ts.cts[0], s.css)
        self._solver = s

    def __dealloc__(self):
        del self.cp


cdef class IC3Bits(__AbstractProver):
    '''
    IC3 variant that splits bit-vectors into individual booleans.
    '''
    def __cinit__(self, Property p, __AbstractTransitionSystem ts, SmtSolver s):
        self.cp = new c_IC3Bits(p.cp[0], ts.cts[0], s.css)
        self._solver = s

    def __dealloc__(self):
        del self.cp


cdef class IC3IA(__AbstractProver):
    '''
    IC3 via Implicit Predicate Abstraction
    '''
    def __cinit__(self, Property p, __AbstractTransitionSystem ts, SmtSolver s):
        self.cp = new c_IC3IA(p.cp[0], ts.cts[0], s.css)
        self._ts = ts
        self._solver = s

    def __dealloc__(self):
        del self.cp


cdef class IC3SA(__AbstractProver):
    '''
    Syntax-Guided Abstraction IC3
    See "Model Checking of Verilog RTL using IC3 with Syntax-Guided Abstraction"
         by Aman Goel and Karem Sakallah in Nasa Formal Methods 2019
    '''
    def __cinit__(self, Property p, __AbstractTransitionSystem ts, SmtSolver s):
        self.cp = new c_IC3SA(p.cp[0], ts.cts[0], s.css)
        self._ts = ts
        self._solver = s

    def __dealloc__(self):
        del self.cp


cdef class InterpolantMC(__AbstractProver):
    def __cinit__(self, Property p, __AbstractTransitionSystem ts, SmtSolver s):
        self.cp = new c_InterpolantMC(p.cp[0], ts.cts[0], s.css)
        self._ts = ts
        self._solver = s

    def __dealloc__(self):
        del self.cp


cdef class ModelBasedIC3(__AbstractProver):
    def __cinit__(self, Property p, __AbstractTransitionSystem ts, SmtSolver s):
        self.cp = new c_ModelBasedIC3(p.cp[0], ts.cts[0], s.css)
        self._ts = ts
        self._solver = s

    def __dealloc__(self):
        del self.cp


IF WITH_MSAT_IC3IA == "ON":
    cdef class MsatIC3IA(__AbstractProver):
        def __cinit__(self, Property p, __AbstractTransitionSystem ts, SmtSolver s):
            self.cp = new c_MsatIC3IA(p.cp[0], ts.cts[0], s.css)
            self._ts = ts
            self._solver = s

        def __dealloc__(self):
            del self.cp


cdef class BTOR2Encoder:
    cdef c_BTOR2Encoder * cbe
    cdef __AbstractTransitionSystem _ts
    def __cinit__(self, str filename, __AbstractTransitionSystem ts):
        self.cbe = new c_BTOR2Encoder(filename.encode(), dref(ts.cts))
        self._ts = ts

    def __dealloc__(self):
        del self.cbe

    def propvec(self):
        res = []
        cdef vector[c_Term] props = dref(self.cbe).propvec()
        for p in props:
            term = Term(self._ts._solver)
            term.ct = p
            res.append(term)
        return res


cdef class SMVEncoder:
    cdef c_SMVEncoder * cbe
    cdef RelationalTransitionSystem _rts
    def __cinit__(self, str filename, RelationalTransitionSystem rts):
        self.cbe = new c_SMVEncoder(filename.encode(), \
                                    dref(<c_RelationalTransitionSystem * ?> rts.cts))
        self._rts = rts

    def __dealloc__(self):
        del self.cbe

    def propvec(self):
        res = []
        cdef vector[c_Term] props = dref(self.cbe).propvec()
        for p in props:
            term = Term(self._rts._solver)
            term.ct = p
            res.append(term)
        return res


cdef class VMTEncoder:
    cdef c_VMTEncoder * cbe
    cdef RelationalTransitionSystem _rts
    def __cinit__(self, str filename, RelationalTransitionSystem rts):
        self.cbe = new c_VMTEncoder(filename.encode(), \
                                    dref(<c_RelationalTransitionSystem * ?> rts.cts))
        self._rts = rts

    def __dealloc__(self):
        del self.cbe

    def propvec(self):
        res = []
        cdef vector[c_Term] props = dref(self.cbe).propvec()
        for p in props:
            term = Term(self._rts._solver)
            term.ct = p
            res.append(term)
        return res


IF WITH_COREIR == "ON":
    cdef class CoreIREncoder:
        cdef c_CoreIREncoder * cbe
        cdef RelationalTransitionSystem _ts
        def __cinit__(self, mod, RelationalTransitionSystem ts):
            self._ts = ts
            cdef uintptr_t adr
            if isinstance(mod, str):
                self.cbe = new c_CoreIREncoder((<string?> (mod.encode())), dref((<c_RelationalTransitionSystem *> ts.cts)))
            elif hasattr(mod, "ptr"):
                adr = <uintptr_t> ctypes.addressof(mod.ptr.contents)
                self.cbe = new c_CoreIREncoder((<c_Module *> adr), dref((<c_RelationalTransitionSystem *> ts.cts)))
            else:
                raise ValueError("CoreIR encoder takes a pycoreir Context or a filename but got {}".format(mod))

        def __dealloc__(self):
            del self.cbe

def add_prop_monitor(__AbstractTransitionSystem ts, Term prop):
    c_add_prop_monitor(dref(ts.cts), prop.ct)

cdef class PonoOptions:
    cdef c_PonoOptions cpo
    def __cinit__(self):
        pass

def parse_options(list opts):
    cdef vector[string] vec
    for o in opts:
        vec.push_back((<str?> o).encode())
    cdef c_PonoOptions c_pono_opts;
    c_pono_opts.parse_and_set_options(vec, False)
    pono_opts = PonoOptions()
    pono_opts.cpo = c_pono_opts

    return pono_opts

cdef class HistoryModifier:
    cdef c_HistoryModifier * chm
    cdef __AbstractTransitionSystem _ts
    def __cinit__(self, __AbstractTransitionSystem ts):
        self.chm = new c_HistoryModifier(dref(ts.cts));
        self._ts = ts

    def __dealloc__(self):
        del self.chm

    def get_hist(self, Term target, int delay):
        if delay < 0:
            raise ValueError("Got negative delay in get_hist: {}".format(delay))

        cdef Term term = Term(self._ts.solver)
        term.ct = dref(self.chm).get_hist(target.ct, delay)
        return term

cdef class VCDWitnessPrinter:
    cdef c_VCDWitnessPrinter * cvwp
    # need to keep the c_cex around on the heap because VCDWitnessPrinter only keeps a reference to it
    cdef vector[c_UnorderedTermMap] * c_cex
    def __cinit__(self, __AbstractTransitionSystem ts, cex):
        c_cex = new vector[c_UnorderedTermMap]()
        for i in range(len(cex)):
            dref(c_cex).push_back(c_UnorderedTermMap())
            for k, v in cex[i].items():
                dref(c_cex)[i][(<Term?> k).ct] = (<Term?> v).ct
        assert len(cex) == dref(c_cex).size(), 'expecting C++ view of witness to be the same length'
        self.cvwp = new c_VCDWitnessPrinter(dref(ts.cts), dref(c_cex))

    def __dealloc__(self):
        del self.cvwp
        del self.c_cex

    def dump_trace_to_file(self, str vcd_file_name):
        dref(self.cvwp).dump_trace_to_file(vcd_file_name.encode())

# standalone functions

def pseudo_init_and_prop(__AbstractTransitionSystem ts, Term prop):
    rts = RelationalTransitionSystem(ts._solver)
    # need to delete old pointer
    del rts.cts
    rts.cts = new c_RelationalTransitionSystem(
        c_pseudo_init_and_prop(dref(ts.cts), prop.ct))
    return rts

def prop_in_trans(__AbstractTransitionSystem ts, Term prop):
    c_prop_in_trans(dref(ts.cts), prop.ct)

def set_global_logger_verbosity(int v):
    c_set_global_logger_verbosity(v)

def check_invar(__AbstractTransitionSystem ts, Term prop, Term invar):
    return c_check_invar(dref(ts.cts), prop.ct, invar.ct)

def coi_reduction(__AbstractTransitionSystem ts, to_keep, verbosity=1):
    '''

    ts - a transition system
    to_keep - a list of terms from ts to keep in the reduction
    Run cone-of-influence reduction on the TransitionSystem ts in-place,
    keeping all variables in the terms in to_keep and any variables that
    influence them in the transition relation.

    '''
    cdef vector[c_Term] c_to_keep
    for t in to_keep:
        c_to_keep.push_back((<Term?> t).ct)

    assert len(to_keep) == c_to_keep.size()
    c_StaticConeOfInfluence(dref(ts.cts), c_to_keep, verbosity)
