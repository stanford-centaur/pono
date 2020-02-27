from libcpp.string cimport string
from libcpp.unordered_map cimport unordered_map
from libcpp.unordered_set cimport unordered_set
from libcpp.vector cimport vector

from smt_switch cimport c_Sort, c_Term, c_SmtSolver, c_TermVec, c_UnorderedTermMap

ctypedef unordered_set[c_Term] c_UnorderedTermSet


cdef extern from "core/ts.h" namespace "cosa":
    cdef cppclass TransitionSystem:
        TransitionSystem(c_SmtSolver & s) except +
        void set_init(const c_Term & init) except +
        void constrain_init(const c_Term & constraint) except +
        void assign_next(const c_Term & state, const c_Term & val) except +
        void add_invar(const c_Term & constraint) except +
        void constrain_inputs(const c_Term & constraint) except +
        void add_constraint(const c_Term & constraint) except +
        void name_term(const string name, const c_Term & t) except +
        c_Term make_input(const string name, const c_Sort & sort) except +
        c_Term make_state(const string name, const c_Sort & sort) except +
        c_Term curr(const c_Term & term) except +
        c_Term next(const c_Term & term) except +
        bint is_curr_var(const c_Term & sv) except +
        bint is_next_var(const c_Term & sv) except +
        c_SmtSolver & solver() except +
        const c_UnorderedTermSet & states() except +
        const c_UnorderedTermSet & inputs() except +
        c_Term init() except +
        c_Term trans() except +
        const c_UnorderedTermMap & state_updates() except +
        unordered_map[string, c_Term] & named_terms() except +
        bint is_functional() except +


cdef extern from "core/rts.h" namespace "cosa":
    cdef cppclass RelationalTransitionSystem(TransitionSystem):
        RelationalTransitionSystem(c_SmtSolver & s) except +
        void set_behavior(const c_Term & init, const c_Term & trans) except +
        void set_trans(const c_Term & trans) except +
        void constrain_trans(const c_Term & constraint) except +


cdef extern from "core/fts.h" namespace "cosa":
    cdef cppclass FunctionalTransitionSystem(TransitionSystem):
        FunctionalTransitionSystem(c_SmtSolver & s) except +


cdef extern from "core/prop.h" namespace "cosa":
    cdef cppclass Property:
        Property(const TransitionSystem& ts, c_Term p) except +
        const c_Term prop() except +
        const TransitionSystem & transition_system() except +


cdef extern from "core/unroller.h" namespace "cosa":
    cdef cppclass Unroller:
        Unroller(const TransitionSystem & ts, c_SmtSolver & solver) except +
        c_Term at_time(const c_Term & t, unsigned int k) except +
        c_Term untime(const c_Term & t) except +


cdef extern from "core/proverresult.h" namespace "cosa":
    cdef cppclass ProverResult:
        pass

    cdef ProverResult UNKNOWN
    cdef ProverResult FALSE
    cdef ProverResult TRUE


cdef extern from "engines/prover.h" namespace "cosa":
    cdef cppclass Prover:
        Prover(const Property & p, c_SmtSolver & s) except +
        void initialize() except +
        ProverResult check_until(int k) except +
        bint witness(vector[c_UnorderedTermMap] & out) except +
        ProverResult prove() except +


cdef extern from "engines/bmc.h" namespace "cosa":
    cdef cppclass Bmc(Prover):
        Bmc(const Property & p, c_SmtSolver & solver) except +


cdef extern from "engines/kinduction.h" namespace "cosa":
    cdef cppclass KInduction(Prover):
        KInduction(const Property & p, c_SmtSolver & solver) except +


cdef extern from "engines/bmc_simplepath.h" namespace "cosa":
    cdef cppclass BmcSimplePath(KInduction):
        BmcSimplePath(const Property & p, c_SmtSolver & solver) except +


cdef extern from "engines/interpolantmc.h" namespace "cosa":
    cdef cppclass InterpolantMC(Prover):
        InterpolantMC(const Property & p, c_SmtSolver & s, c_SmtSolver & interpolator) except +


cdef extern from "frontends/btor2_encoder.h" namespace "cosa":
    cdef cppclass BTOR2Encoder:
        pass
