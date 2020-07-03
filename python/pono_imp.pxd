from libcpp.string cimport string
from libcpp.unordered_map cimport unordered_map
from libcpp.unordered_set cimport unordered_set
from libcpp.vector cimport vector

from smt_switch cimport c_Sort, c_Term, c_SmtSolver, c_TermVec, c_UnorderedTermMap

ctypedef unordered_set[c_Term] c_UnorderedTermSet


cdef extern from "core/ts.h" namespace "pono":
    cdef cppclass TransitionSystem:
        TransitionSystem(c_SmtSolver & s) except +
        void set_init(const c_Term & init) except +
        void constrain_init(const c_Term & constraint) except +
        void assign_next(const c_Term & state, const c_Term & val) except +
        void add_invar(const c_Term & constraint) except +
        void constrain_inputs(const c_Term & constraint) except +
        void add_constraint(const c_Term & constraint) except +
        void name_term(const string name, const c_Term & t) except +
        c_Term make_inputvar(const string name, const c_Sort & sort) except +
        c_Term make_statevar(const string name, const c_Sort & sort) except +
        c_Term curr(const c_Term & term) except +
        c_Term next(const c_Term & term) except +
        bint is_curr_var(const c_Term & sv) except +
        bint is_next_var(const c_Term & sv) except +
        c_SmtSolver & solver() except +
        const c_UnorderedTermSet & statevars() except +
        const c_UnorderedTermSet & inputvars() except +
        c_Term init() except +
        c_Term trans() except +
        const c_UnorderedTermMap & state_updates() except +
        unordered_map[string, c_Term] & named_terms() except +
        bint is_functional() except +


cdef extern from "core/rts.h" namespace "pono":
    cdef cppclass RelationalTransitionSystem(TransitionSystem):
        RelationalTransitionSystem(c_SmtSolver & s) except +
        void set_behavior(const c_Term & init, const c_Term & trans) except +
        void set_trans(const c_Term & trans) except +
        void constrain_trans(const c_Term & constraint) except +


cdef extern from "core/fts.h" namespace "pono":
    cdef cppclass FunctionalTransitionSystem(TransitionSystem):
        FunctionalTransitionSystem(c_SmtSolver & s) except +


cdef extern from "core/prop.h" namespace "pono":
    cdef cppclass Property:
        Property(const TransitionSystem& ts, c_Term p) except +
        const c_Term prop() except +
        const TransitionSystem & transition_system() except +


cdef extern from "core/unroller.h" namespace "pono":
    cdef cppclass Unroller:
        Unroller(const TransitionSystem & ts, c_SmtSolver & solver) except +
        c_Term at_time(const c_Term & t, unsigned int k) except +
        c_Term untime(const c_Term & t) except +


cdef extern from "core/proverresult.h" namespace "pono":
    cdef cppclass ProverResult:
        pass

    cdef ProverResult UNKNOWN
    cdef ProverResult FALSE
    cdef ProverResult TRUE


cdef extern from "engines/prover.h" namespace "pono":
    cdef cppclass Prover:
        Prover(const Property & p, c_SmtSolver & s) except +
        void initialize() except +
        ProverResult check_until(int k) except +
        bint witness(vector[c_UnorderedTermMap] & out) except +
        ProverResult prove() except +


cdef extern from "engines/bmc.h" namespace "pono":
    cdef cppclass Bmc(Prover):
        Bmc(const Property & p, c_SmtSolver & solver) except +


cdef extern from "engines/kinduction.h" namespace "pono":
    cdef cppclass KInduction(Prover):
        KInduction(const Property & p, c_SmtSolver & solver) except +


cdef extern from "engines/bmc_simplepath.h" namespace "pono":
    cdef cppclass BmcSimplePath(KInduction):
        BmcSimplePath(const Property & p, c_SmtSolver & solver) except +


cdef extern from "engines/interpolantmc.h" namespace "pono":
    cdef cppclass InterpolantMC(Prover):
        InterpolantMC(const Property & p, c_SmtSolver & s, c_SmtSolver & interpolator) except +


cdef extern from "frontends/btor2_encoder.h" namespace "pono":
    cdef cppclass BTOR2Encoder:
        BTOR2Encoder(string filename, TransitionSystem & ts) except +


# WITH_COREIR is set in python/CMakeLists.txt via the --compile-time-env flag of Cython
IF WITH_COREIR == "ON":
    cdef extern from "frontends/coreir_encoder.h" namespace "pono":
        cdef cppclass CoreIREncoder:
            CoreIREncoder(string filename, RelationalTransitionSystem & ts) except +


cdef extern from "utils/logger.h" namespace "pono":
    void set_global_logger_verbosity(unsigned int v) except +


cdef extern from "utils/term_analysis.h" namespace "pono":
    void get_free_symbols(const c_Term & term, c_UnorderedTermSet & out_symbols) except +
