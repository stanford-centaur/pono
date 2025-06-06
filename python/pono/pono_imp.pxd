from libc.stdint cimport uint64_t
from libcpp cimport bool
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.unordered_map cimport unordered_map
from libcpp.unordered_set cimport unordered_set
from libcpp.vector cimport vector

from smt_switch cimport c_SortKind, c_Sort, c_PrimOp, c_Op, \
    c_Term, c_SmtSolver, c_SortVec, c_TermVec, c_UnorderedTermMap

ctypedef unordered_set[c_Term] c_UnorderedTermSet


cdef extern from "core/ts.h" namespace "pono":
    cdef cppclass TransitionSystem:
        TransitionSystem() except +
        TransitionSystem(c_SmtSolver & s) except +
        void set_init(const c_Term & init) except +
        void constrain_init(const c_Term & constraint) except +
        void assign_next(const c_Term & state, const c_Term & val) except +
        void add_invar(const c_Term & constraint) except +
        void constrain_inputs(const c_Term & constraint) except +
        void add_constraint(const c_Term & constraint, bint to_init_and_next) except +
        void name_term(const string name, const c_Term & t) except +
        c_Term make_inputvar(const string name, const c_Sort & sort) except +
        c_Term make_statevar(const string name, const c_Sort & sort) except +
        c_Term curr(const c_Term & term) except +
        c_Term next(const c_Term & term) except +
        bint is_curr_var(const c_Term & sv) except +
        bint is_next_var(const c_Term & sv) except +
        bint is_input_var(const c_Term & sv) except +
        string get_name(const c_Term & t) except +
        c_Term lookup(string name) except +
        void add_statevar(const c_Term & cv, const c_Term & nv) except +
        void add_inputvar(const c_Term & v) except +
        c_SmtSolver & solver() except +
        const c_UnorderedTermSet & statevars() except +
        const c_UnorderedTermSet & inputvars() except +
        c_Term init() except +
        c_Term trans() except +
        const c_UnorderedTermMap & state_updates() except +
        unordered_map[string, c_Term] & named_terms() except +
        const vector[pair[c_Term, bool]] & constraints() except +
        bint is_functional() except +
        bint is_deterministic() except +
        void drop_state_updates(const c_TermVec & svs) except +
        void promote_inputvar(const c_Term & iv) except +
        void replace_terms(const c_UnorderedTermMap & to_replace) except +
        c_Sort make_sort(const string name, uint64_t arity) except +
        c_Sort make_sort(const c_SortKind sk) except +
        c_Sort make_sort(const c_SortKind sk, uint64_t size) except +
        c_Sort make_sort(const c_SortKind sk, const c_Sort & sort1) except +
        c_Sort make_sort(const c_SortKind sk, const c_Sort & sort1, const c_Sort & sort2) except +
        c_Sort make_sort(const c_SortKind sk, const c_Sort & sort1, const c_Sort & sort2, const c_Sort & sort3) except +
        c_Sort make_sort(const c_SortKind sk, const c_SortVec & sorts) except +
        c_Term make_term(bint b) except +
        c_Term make_term(const string val, const c_Sort & sort) except +
        c_Term make_term(const string val, const c_Sort & sort, uint64_t base) except +
        c_Term make_term(const c_Term & val, const c_Sort & sort) except +
        c_Term make_term(const c_Op op, const c_TermVec & terms) except +


cdef extern from "core/rts.h" namespace "pono":
    cdef cppclass RelationalTransitionSystem(TransitionSystem):
        RelationalTransitionSystem(c_SmtSolver & s) except +
        RelationalTransitionSystem(const TransitionSystem & rts) except +
        void set_behavior(const c_Term & init, const c_Term & trans) except +
        void set_trans(const c_Term & trans) except +
        void constrain_trans(const c_Term & constraint) except +


cdef extern from "core/fts.h" namespace "pono":
    cdef cppclass FunctionalTransitionSystem(TransitionSystem):
        FunctionalTransitionSystem(c_SmtSolver & s) except +
        FunctionalTransitionSystem(const TransitionSystem & fts) except +


cdef extern from "core/prop.h" namespace "pono":
    cdef cppclass Property:
        Property(const c_SmtSolver& s, c_Term p) except +
        const c_Term prop() except +
        const c_SmtSolver & solver() except +
        string name()

cdef extern from "core/unroller.h" namespace "pono":
    cdef cppclass Unroller:
        Unroller(const TransitionSystem & ts, const string & time_id) except +
        c_Term at_time(const c_Term & t, unsigned int k) except +
        c_Term untime(const c_Term & t) except +
        size_t get_var_time(const c_Term & v) except +


cdef extern from "core/proverresult.h" namespace "pono":
    cdef cppclass ProverResult:
        pass

    cdef ProverResult UNKNOWN
    cdef ProverResult FALSE
    cdef ProverResult TRUE


cdef extern from "engines/prover.h" namespace "pono":
    cdef cppclass Prover:
        Prover(const Property & p, const TransitionSystem & ts,
               c_SmtSolver & s) except +
        void initialize() except +
        ProverResult check_until(int k) except +
        bint witness(vector[c_UnorderedTermMap] & out) except +
        c_Term invar() except +
        ProverResult prove() except +


cdef extern from "engines/bmc.h" namespace "pono":
    cdef cppclass Bmc(Prover):
        Bmc(const Property & p, const TransitionSystem & ts,
            c_SmtSolver & solver) except +


cdef extern from "engines/kinduction.h" namespace "pono":
    cdef cppclass KInduction(Prover):
        KInduction(const Property & p, const TransitionSystem & ts,
                   c_SmtSolver & solver) except +


cdef extern from "engines/bmc_simplepath.h" namespace "pono":
    cdef cppclass BmcSimplePath(KInduction):
        BmcSimplePath(const Property & p, const TransitionSystem & ts,
                      c_SmtSolver & solver) except +


cdef extern from "engines/ic3.h" namespace "pono":
    cdef cppclass IC3(Prover):
        IC3(const Property & p, const TransitionSystem & ts,
            c_SmtSolver & solver) except +


cdef extern from "engines/ic3bits.h" namespace "pono":
    cdef cppclass IC3Bits(Prover):
        IC3Bits(const Property & p, const TransitionSystem & ts,
                c_SmtSolver & solver) except +


cdef extern from "engines/ic3ia.h" namespace "pono":
    cdef cppclass IC3IA(Prover):
        IC3IA(const Property & p, const TransitionSystem & ts,
              c_SmtSolver & solver) except +


cdef extern from "engines/ic3sa.h" namespace "pono":
    cdef cppclass IC3SA(Prover):
        IC3SA(const Property & p, const TransitionSystem & ts,
              c_SmtSolver & solver) except +


cdef extern from "engines/interpolantmc.h" namespace "pono":
    cdef cppclass InterpolantMC(Prover):
        InterpolantMC(const Property & p, const TransitionSystem & ts,
                      c_SmtSolver & s) except +


cdef extern from "engines/mbic3.h" namespace "pono":
    cdef cppclass ModelBasedIC3(Prover):
        ModelBasedIC3(const Property & p, const TransitionSystem & ts,
                      c_SmtSolver & solver) except +


# WITH_MSAT_IC3IA is set in python/CMakeLists.txt via the --compile-time-env flag of Cython
IF WITH_MSAT_IC3IA == "ON":
    cdef extern from "engines/msat_ic3ia.h" namespace "pono":
        cdef cppclass MsatIC3IA(Prover):
            MsatIC3IA(const Property & p, const TransitionSystem & ts,
                      c_SmtSolver & s) except +


cdef extern from "frontends/btor2_encoder.h" namespace "pono":
    cdef cppclass BTOR2Encoder:
        BTOR2Encoder(string filename, TransitionSystem & ts) except +
        const c_TermVec & propvec() except +


cdef extern from "frontends/smv_encoder.h" namespace "pono":
    cdef cppclass SMVEncoder:
        SMVEncoder(string filename, RelationalTransitionSystem & ts) except +
        const c_TermVec & propvec() except +


cdef extern from "frontends/vmt_encoder.h" namespace "pono":
    cdef cppclass VMTEncoder:
        VMTEncoder(string filename, RelationalTransitionSystem & ts) except +
        const c_TermVec & propvec() except +


# WITH_COREIR is set in python/CMakeLists.txt via the --compile-time-env flag of Cython
IF WITH_COREIR == "ON":
    cdef extern from "coreir.h" namespace "CoreIR":
        cdef cppclass Module:
            pass

    cdef extern from "frontends/coreir_encoder.h" namespace "pono":
        cdef cppclass CoreIREncoder:
            CoreIREncoder(string filename, RelationalTransitionSystem & ts) except +
            CoreIREncoder(Module * top_mod, RelationalTransitionSystem & ts) except +


cdef extern from "modifiers/static_coi.h" namespace "pono":
    cdef cppclass StaticConeOfInfluence:
        StaticConeOfInfluence(TransitionSystem & ts,
                        const c_TermVec & to_keep,
                        int verbosity) except +

cdef extern from "modifiers/prop_monitor.h" namespace "pono":
    c_Term add_prop_monitor(TransitionSystem & ts,
                            const c_Term & prop) except +

cdef extern from "modifiers/history_modifier.h" namespace "pono":
    cdef cppclass HistoryModifier:
        HistoryModifier(TransitionSystem & ts) except +
        c_Term get_hist(const c_Term & target, size_t delay) except +

cdef extern from "modifiers/mod_ts_prop.h" namespace "pono":
    TransitionSystem pseudo_init_and_prop(TransitionSystem & ts, c_Term & prop) except +
    void prop_in_trans(TransitionSystem & ts, const c_Term & prop) except +

cdef extern from "options/options.h" namespace "pono":
    cdef cppclass PonoOptions:
        PonoOptions() except +
        ProverResult parse_and_set_options(vector[string] & opts, bint expect_file) except +

cdef extern from "printers/vcd_witness_printer.h" namespace "pono":
    cdef cppclass VCDWitnessPrinter:
        VCDWitnessPrinter(const TransitionSystem & ts,
                          vector[c_UnorderedTermMap] & cex)
        void dump_trace_to_file(const string & vcd_file_name) except +

cdef extern from "utils/logger.h" namespace "pono":
    void set_global_logger_verbosity(unsigned int v) except +

cdef extern from "utils/ts_analysis.h" namespace "pono":
    bint check_invar(const TransitionSystem & ts,
                     const c_Term & prop,
                     const c_Term & invar) except +

