cmake_minimum_required(VERSION 3.18...4.0)

project(pono VERSION 2.1.0)

# Pre-release/dev suffix appended to PROJECT_VERSION. project(VERSION) accepts
# only numeric components, so prerelease markers live here. Spelled in PEP 440
# form ("b1", "rc1", ".dev0") so that the C++ and Python versions stay a single
# string. Empty for final releases.
set(PONO_VERSION_SUFFIX ".dev0")

# The declared version, and the single source of truth for every version
# pono reports: --version, the wheel, and pono.__version__.
set(PONO_RELEASE_VERSION "${PROJECT_VERSION}${PONO_VERSION_SUFFIX}")

# Generate compilation database compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include("${PROJECT_SOURCE_DIR}/cmake/PonoVersion.cmake")

# handle different versions of CMake
if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.0 AND NOT APPLE)
  set(CMAKE_CXX17_STANDARD_COMPILE_OPTION "-std=c++17")
  set(CMAKE_CXX17_EXTENSION_COMPILE_OPTION "-std=gnu++17")
elseif(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.1 OR APPLE)
  set(CMAKE_CXX17_STANDARD_COMPILE_OPTION "-std=c++1z")
  set(CMAKE_CXX17_EXTENSION_COMPILE_OPTION "-std=gnu++1z")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_definitions(-DPONO_SRC_DIR=${PROJECT_SOURCE_DIR})

# Whether pono-bin is actually linked with -static. PONO_STATIC_EXEC asks for
# that, but static libraries cannot be built on Macs, so the request only
# draws a warning there. The fmt provisioning below and the link flag on
# pono-bin both read this instead of repeating the platform test.
# See: https://developer.apple.com/library/archive/qa/qa1118/_index.html
set(PONO_LINK_STATIC OFF)

if(PONO_STATIC_EXEC STREQUAL "YES")
  if(WITH_COREIR)
    message(
      FATAL_ERROR
      "CoreIR support and a fully static executable are incompatible, since "
      "CoreIR is only available as a shared library. Disable one or the other."
    )
  endif()
  if(NOT PONO_LIB_TYPE STREQUAL "STATIC")
    message(
      VERBOSE
      "PONO_STATIC_EXEC=YES requires a static pono-lib; overriding "
      "PONO_LIB_TYPE=${PONO_LIB_TYPE} to STATIC."
    )
    set(PONO_LIB_TYPE STATIC)
  endif()
  if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    message(WARNING "Static executables are not supported on macOS")
  else()
    set(PONO_LINK_STATIC ON)
  endif()
endif()

# Defines pono-warnings, the interface target holding the compiler warning
# flags, which pono's own targets link against below.
include("${PROJECT_SOURCE_DIR}/cmake/PonoWarnings.cmake")

if(WITH_BOOLECTOR)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_BOOLECTOR")
endif()

if(WITH_MSAT)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_MSAT")
endif()

if(WITH_YICES2)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_YICES2")
endif()

if(WITH_Z3)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_Z3")
endif()

if(WITH_PROFILING)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_PROFILING")
endif()

# Try to find our custom bison build, falling back to Homebrew.
if(IS_DIRECTORY "${PROJECT_SOURCE_DIR}/deps/bison-install")
  set(BISON_ROOT "${PROJECT_SOURCE_DIR}/deps/bison-install")
else()
  execute_process(
    COMMAND brew --prefix
    RESULT_VARIABLE HOMEBREW_RESULT
    OUTPUT_VARIABLE HOMEBREW_PREFIX
    OUTPUT_QUIET
    ERROR_QUIET
  )
  if(HOMEBREW_RESULT EQUAL 0 AND IS_DIRECTORY "${HOMEBREW_PREFIX}/opt/bison")
    set(BISON_ROOT "${HOMEBREW_PREFIX}/opt/bison")
  endif()
endif()
find_package(BISON 3.4.2 REQUIRED)
find_package(FLEX 2.6.4 REQUIRED)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# Install prefixes of the dependencies pono locates itself, one per Find
# module in cmake/. Declared with an empty default so that configure.sh stays
# the one place defaults live, and so that leaving a backend disabled does not
# warn about an unused -D setting. IC3IA_ROOT is a source tree rather than an
# install prefix, since ic3ia has no install step.
set(Bitwuzla_ROOT "" CACHE PATH "bitwuzla install prefix")
set(Btor2Tools_ROOT "" CACHE PATH "Btor2Tools install prefix")
set(CoreIR_ROOT "" CACHE PATH "CoreIR install prefix")
set(IC3IA_ROOT "" CACHE PATH "ic3ia source tree")
set(MathSAT_ROOT "" CACHE PATH "MathSAT install prefix")
set(SmtSwitch_ROOT "" CACHE PATH "smt-switch install prefix")
set(Z3_ROOT "" CACHE PATH "Z3 install prefix")

find_package(GMPXX REQUIRED)
find_package(Btor2Tools REQUIRED)

# Defines fmt::fmt-header-only, from a system fmt or a fetched one.
include("${PROJECT_SOURCE_DIR}/cmake/PonoFmt.cmake")

if(WITH_PROFILING)
  find_package(Gperftools REQUIRED)
endif()

if(WITH_COREIR)
  find_package(CoreIR REQUIRED)
endif()

if(WITH_MSAT OR WITH_MSAT_IC3IA)
  # pono needs direct access to MathSAT for setting solver options.
  find_package(MathSAT REQUIRED)
endif()

if(WITH_MSAT_IC3IA)
  find_package(IC3IA REQUIRED)
endif()

set(SmtSwitch_FIND_COMPONENTS bitwuzla cvc5)
if(WITH_BOOLECTOR)
  list(APPEND SmtSwitch_FIND_COMPONENTS btor)
endif()
if(WITH_MSAT)
  list(APPEND SmtSwitch_FIND_COMPONENTS msat)
endif()
if(WITH_YICES2)
  list(APPEND SmtSwitch_FIND_COMPONENTS yices2)
endif()
if(WITH_Z3)
  list(APPEND SmtSwitch_FIND_COMPONENTS z3)
endif()
find_package(SmtSwitch REQUIRED COMPONENTS ${SmtSwitch_FIND_COMPONENTS})

# Generate the SMV parser using Flex/Bison
bison_target(
  SMVParser
  ${PROJECT_SOURCE_DIR}/frontends/smvparser.y
  ${CMAKE_CURRENT_BINARY_DIR}/smvparser.cpp
  DEFINES_FILE ${CMAKE_CURRENT_BINARY_DIR}/smvparser.h
)

flex_target(
  SMVScanner
  ${PROJECT_SOURCE_DIR}/frontends/smvparser.l
  ${CMAKE_CURRENT_BINARY_DIR}/smvlexer.cpp
)

# Warnings the Flex and Bison skeletons raise about their own output, which
# the rules in smvparser.l and smvparser.y cannot influence. The positive
# forms are what get probed: a compiler without the warning reports the -Wno-
# as unrecognized rather than ignoring it.
#
# The parser's unreachable breaks are the one case the grammar has a hand in:
# Bison ends every action with a break, and the actions for the constructs
# pono rejects end in a throw. Reaching the break would mean accepting them.
pono_filter_supported_flags(
  lexer_warnings
  -Wextra-semi-stmt # a ; that forms an empty statement
  -Wmisleading-indentation # a statement indented under an if it is not in
)
pono_filter_supported_flags(
  parser_warnings
  -Wunreachable-code-break # a break after an action that throws
  -Wunused-but-set-variable # yynerrs_, which the skeleton never reads
)
list(TRANSFORM lexer_warnings REPLACE "^-W" "-Wno-")
list(TRANSFORM parser_warnings REPLACE "^-W" "-Wno-")
set_source_files_properties(
  "${CMAKE_CURRENT_BINARY_DIR}/smvlexer.cpp"
  PROPERTIES COMPILE_OPTIONS "${lexer_warnings}"
)
set_source_files_properties(
  "${CMAKE_CURRENT_BINARY_DIR}/smvparser.cpp"
  PROPERTIES COMPILE_OPTIONS "${parser_warnings}"
)

add_flex_bison_dependency(SMVScanner SMVParser)

set(
  SOURCES
  "${PROJECT_SOURCE_DIR}/core/ts.cpp"
  "${PROJECT_SOURCE_DIR}/core/rts.cpp"
  "${PROJECT_SOURCE_DIR}/core/fts.cpp"
  "${PROJECT_SOURCE_DIR}/core/unroller.cpp"
  "${PROJECT_SOURCE_DIR}/core/functional_unroller.cpp"
  "${PROJECT_SOURCE_DIR}/core/prop.cpp"
  "${PROJECT_SOURCE_DIR}/core/proverresult.cpp"
  "${PROJECT_SOURCE_DIR}/engines/prover.cpp"
  "${PROJECT_SOURCE_DIR}/engines/bmc.cpp"
  "${PROJECT_SOURCE_DIR}/engines/bmc_simplepath.cpp"
  "${PROJECT_SOURCE_DIR}/engines/cegar_ops_uf.cpp"
  "${PROJECT_SOURCE_DIR}/engines/cegar_values.cpp"
  "${PROJECT_SOURCE_DIR}/engines/ceg_prophecy_arrays.cpp"
  "${PROJECT_SOURCE_DIR}/engines/dual_approx_reach.cpp"
  "${PROJECT_SOURCE_DIR}/engines/ic3.cpp"
  "${PROJECT_SOURCE_DIR}/engines/ic3base.cpp"
  "${PROJECT_SOURCE_DIR}/engines/ic3bits.cpp"
  "${PROJECT_SOURCE_DIR}/engines/ic3ia.cpp"
  "${PROJECT_SOURCE_DIR}/engines/ic3sa.cpp"
  "${PROJECT_SOURCE_DIR}/engines/interpolantmc.cpp"
  "${PROJECT_SOURCE_DIR}/engines/interp_seq_mc.cpp"
  "${PROJECT_SOURCE_DIR}/engines/kinduction.cpp"
  "${PROJECT_SOURCE_DIR}/engines/kliveness.cpp"
  "${PROJECT_SOURCE_DIR}/engines/mbic3.cpp"
  "${PROJECT_SOURCE_DIR}/engines/syguspdr.cpp"
  "${PROJECT_SOURCE_DIR}/frontends/btor2_encoder.cpp"
  "${PROJECT_SOURCE_DIR}/frontends/smv_encoder.cpp"
  "${PROJECT_SOURCE_DIR}/frontends/smv_node.cpp"
  "${PROJECT_SOURCE_DIR}/frontends/vmt_encoder.cpp"
  "${PROJECT_SOURCE_DIR}/modifiers/array_abstractor.cpp"
  "${PROJECT_SOURCE_DIR}/modifiers/control_signals.cpp"
  "${PROJECT_SOURCE_DIR}/modifiers/implicit_predicate_abstractor.cpp"
  "${PROJECT_SOURCE_DIR}/modifiers/liveness_to_safety_translator.cpp"
  "${PROJECT_SOURCE_DIR}/modifiers/history_modifier.cpp"
  "${PROJECT_SOURCE_DIR}/modifiers/mod_ts_prop.cpp"
  "${PROJECT_SOURCE_DIR}/modifiers/ops_abstractor.cpp"
  "${PROJECT_SOURCE_DIR}/modifiers/prop_monitor.cpp"
  "${PROJECT_SOURCE_DIR}/modifiers/prophecy_modifier.cpp"
  "${PROJECT_SOURCE_DIR}/modifiers/static_coi.cpp"
  "${PROJECT_SOURCE_DIR}/modifiers/op_abstractor.cpp"
  "${PROJECT_SOURCE_DIR}/printers/btor2_witness_printer.cpp"
  "${PROJECT_SOURCE_DIR}/printers/vcd_witness_printer.cpp"
  "${PROJECT_SOURCE_DIR}/refiners/array_axiom_enumerator.cpp"
  "${PROJECT_SOURCE_DIR}/smt/available_solvers.cpp"
  "${PROJECT_SOURCE_DIR}/utils/fcoi.cpp"
  "${PROJECT_SOURCE_DIR}/utils/logger.cpp"
  "${PROJECT_SOURCE_DIR}/utils/make_provers.cpp"
  "${PROJECT_SOURCE_DIR}/utils/term_analysis.cpp"
  "${PROJECT_SOURCE_DIR}/utils/term_walkers.cpp"
  "${PROJECT_SOURCE_DIR}/utils/ts_analysis.cpp"
  "${PROJECT_SOURCE_DIR}/utils/ts_manipulation.cpp"
  "${PROJECT_SOURCE_DIR}/utils/sygus_ic3formula_helper.cpp"
  "${PROJECT_SOURCE_DIR}/utils/sygus_predicate_constructor.cpp"
  "${PROJECT_SOURCE_DIR}/utils/str_util.cpp"
  "${PROJECT_SOURCE_DIR}/utils/partial_model.cpp"
  "${PROJECT_SOURCE_DIR}/utils/syntax_analysis_common.cpp"
  "${PROJECT_SOURCE_DIR}/utils/syntax_analysis_walker.cpp"
  "${PROJECT_SOURCE_DIR}/utils/syntax_analysis.cpp"
  "${PROJECT_SOURCE_DIR}/options/options.cpp"
  "${BISON_SMVParser_OUTPUTS}"
  "${FLEX_SMVScanner_OUTPUTS}"
)

if(WITH_COREIR)
  add_definitions(-DWITH_COREIR)
  set(SOURCES "${SOURCES}" "${PROJECT_SOURCE_DIR}/frontends/coreir_encoder.cpp")
endif()

if(WITH_MSAT_IC3IA)
  add_definitions(-DWITH_MSAT_IC3IA)
  set(SOURCES "${SOURCES}" "${PROJECT_SOURCE_DIR}/engines/msat_ic3ia.cpp")
endif()

add_library(pono-lib "${PONO_LIB_TYPE}" ${SOURCES})
set_target_properties(pono-lib PROPERTIES OUTPUT_NAME pono)

# A target-level edge, so that both the Ninja and Makefile generators finish
# writing version.h before starting any pono-lib object. Everything else links
# pono-lib and inherits the ordering.
add_dependencies(pono-lib pono-version-header)

target_include_directories(pono-lib PUBLIC "${PROJECT_SOURCE_DIR}")

# Generated and vendored headers, whose warnings are not ours to act on.
# SYSTEM keeps the compiler quiet about them. Consumers pick this up through
# pono-lib. Everything in the build directory is generated: smvparser.h,
# location.hh and version.h.
#
# A header that a system header includes counts as one too, and smvparser.h
# includes frontends/smv_node.h. Warnings in smv_node.h therefore reach us
# only from the sources that include it themselves, which today is enough to
# report every one of them.
target_include_directories(
  pono-lib
  SYSTEM
  PUBLIC "${CMAKE_CURRENT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/contrib/optionparser-1.7/src"
)

target_link_libraries(pono-lib PRIVATE pono-warnings)

# Should we build python bindings
option(BUILD_PYTHON_BINDINGS "Build Python bindings")

if(BUILD_PYTHON_BINDINGS)
  add_subdirectory(python)
endif()

target_link_libraries(pono-lib PUBLIC smt-switch-bitwuzla)
target_link_libraries(pono-lib PUBLIC smt-switch-cvc5)

if(WITH_BOOLECTOR)
  target_link_libraries(pono-lib PUBLIC smt-switch-btor)
endif()

if(WITH_MSAT)
  target_link_libraries(pono-lib PUBLIC smt-switch-msat)
endif()

if(WITH_YICES2)
  target_link_libraries(pono-lib PUBLIC smt-switch-yices2)
endif()

if(WITH_Z3)
  target_link_libraries(pono-lib PRIVATE smt-switch-z3)
endif()

# MathSAT::mathsat carries only an include directory; pono links MathSAT
# itself through smt-switch. Both backends need those headers.
if(WITH_MSAT OR WITH_MSAT_IC3IA)
  target_link_libraries(pono-lib PUBLIC MathSAT::mathsat)
endif()

if(WITH_MSAT_IC3IA)
  target_link_libraries(pono-lib PUBLIC IC3IA::ic3ia)
endif()

if(WITH_COREIR)
  target_link_libraries(pono-lib PUBLIC CoreIR::coreir CoreIR::verilogAST)
endif()

target_link_libraries(pono-lib PUBLIC smt-switch)
target_link_libraries(pono-lib PUBLIC Btor2Tools::btor2parser)
target_link_libraries(pono-lib PUBLIC GMPXX::gmpxx)
# PUBLIC because the installed utils/logger.h includes <fmt/format.h>. The
# header-only target carries fmt's include directory without a library to link,
# so this stays compatible with PONO_STATIC_EXEC.
target_link_libraries(pono-lib PUBLIC fmt::fmt-header-only)

if(WITH_PROFILING)
  target_link_libraries(pono-lib PUBLIC Gperftools::profiler)
endif()

add_executable(pono-bin "${PROJECT_SOURCE_DIR}/pono.cpp")
set_target_properties(pono-bin PROPERTIES OUTPUT_NAME pono)

target_include_directories(
  pono-bin
  PUBLIC "${PROJECT_SOURCE_DIR}/utils" "${PROJECT_SOURCE_DIR}/core"
)

target_link_libraries(pono-bin PUBLIC pono-lib)
target_link_libraries(pono-bin PRIVATE pono-warnings)

if(PONO_LINK_STATIC)
  target_link_libraries(pono-bin PUBLIC -static)
endif()

enable_testing()
# Add tests subdirectory
# The CMakeLists.txt file there sets up googletest
# and builds all the parametrized tests
add_subdirectory(tests EXCLUDE_FROM_ALL)

# Documentation
if(BUILD_DOCS)
  find_package(Doxygen)
  doxygen_add_docs(
    docs
    ALL
    core/
    engines/
    frontends/
    modifiers/
    options/
    printers/
    refiners/
    smt/
    utils/
    pono.cpp
  )
endif()

set(CMAKE_INSTALL_MESSAGE LAZY)

# install smt-switch
install(TARGETS pono-lib DESTINATION lib)
install(TARGETS pono-bin DESTINATION bin)

# install public headers
install(
  DIRECTORY "${PROJECT_SOURCE_DIR}/core/"
  DESTINATION include/pono/core
  FILES_MATCHING
  PATTERN "*.h"
)

install(
  DIRECTORY "${PROJECT_SOURCE_DIR}/engines/"
  DESTINATION include/pono/engines
  FILES_MATCHING
  PATTERN "*.h"
)

install(
  DIRECTORY "${PROJECT_SOURCE_DIR}/frontends/"
  DESTINATION include/pono/frontends
  FILES_MATCHING
  PATTERN "*.h"
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/location.hh" DESTINATION include/pono/)

install(
  DIRECTORY "${PROJECT_SOURCE_DIR}/modifiers/"
  DESTINATION include/pono/modifiers
  FILES_MATCHING
  PATTERN "*.h"
)

install(
  DIRECTORY "${PROJECT_SOURCE_DIR}/options/"
  DESTINATION include/pono/options
  FILES_MATCHING
  PATTERN "*.h"
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/smvparser.h" DESTINATION include/pono/)

install(
  DIRECTORY "${PROJECT_SOURCE_DIR}/utils/"
  DESTINATION include/pono/utils
  FILES_MATCHING
  PATTERN "*.h"
)

# uninstall target
# copied from https://gitlab.kitware.com/cmake/community/wikis/FAQ#can-i-do-make-uninstall-with-cmake
if(NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE
    @ONLY
  )

  add_custom_target(
    uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
  )
endif()
