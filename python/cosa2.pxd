from libcpp.string cimport string
from libcpp.unordered_map cimport unordered_map

from smt_switch_imp cimport Sort, Term, SmtSolver, TermVec, UnorderedTermMap, UnorderedTermSet


cdef extern from "ts.h" namespace "cosa":
    cdef cppclass TransitionSystem:
        TransitionSystem(SmtSolver & s) except +
        void set_behavior(const Term & init, const Term & trans) except +
        void set_init(const Term & init) except +
        void constrain_init(const Term & constraint) except +
        void set_trans(const Term & trans) except +
        void constrain_trans(const Term & constraint) except +
        void assign_next(const Term & state, const Term & val) except +
        void add_invar(const Term & constraint) except +
        void constrain_inputs(const Term & constraint) except +
        void add_constraint(const Term & constraint) except +
        void name_term(const string name, const Term & t) except +
        Term make_input(const string name, const Sort & sort) except +
        Term make_state(const string name, const Sort & sort) except +
        Term curr(const Term & term) except +
        Term next(const Term & term) except +
        Term is_curr_var(const Term & sv) except +
        Term is_next_var(const Term & sv) except +
        SmtSolver & solver() except +
        const UnorderedTermSet & states() except +
        const UnorderedTermSet & inputs() except +
        Term init() except +
        Term trans() except +
        const UnorderedTermMap & state_updates() except +
        # TODO: add this back in
        # unordered_map<string, Term> & named_terms() except +
        bint is_functional() except +


cdef extern from "fts.h" namespace "cosa":
    cdef cppclass FunctionalTransitionSystem:
        pass


cdef extern from "prop.h" namespace "cosa":
    cdef cppclass Property:
        pass


cdef extern from "unroller.h" namespace "cosa":
    cdef cppclass Unroller:
        pass


cdef extern from "bmc.h" namespace "cosa":
    cdef cppclass Bmc:
        pass


cdef extern from "kinduction.h" namespace "cosa":
    cdef cppclass KInduction:
        pass


cdef extern from "interpolantmc.h" namespace "cosa":
    cdef cppclass InterpolantMC:
        pass
