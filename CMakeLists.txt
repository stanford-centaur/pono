cmake_minimum_required(VERSION 3.1)

project(pono)

set(PONO_MAJOR 0)   # Major component of Pono version
set(PONO_MINOR 1)   # Minor component of Pono version
set(PONO_RELEASE 1) # Release component of Pono version

# handle different versions of CMake
if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.0 AND NOT APPLE)
  set(CMAKE_CXX17_STANDARD_COMPILE_OPTION "-std=c++17")
  set(CMAKE_CXX17_EXTENSION_COMPILE_OPTION "-std=gnu++17")
elseif (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.1 OR APPLE)
  set(CMAKE_CXX17_STANDARD_COMPILE_OPTION "-std=c++1z")
  set(CMAKE_CXX17_EXTENSION_COMPILE_OPTION "-std=gnu++1z")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_definitions(-DPONO_SRC_DIR=${PROJECT_SOURCE_DIR})

if(APPLE)
  set(SHARED_LIB_EXT "dylib")
else()
  set(SHARED_LIB_EXT "so")
endif()

if (WITH_MSAT)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_MSAT")
endif()

if (WITH_PROFILING)
  find_library(GOOGLE_PERF profiler REQUIRED)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_PROFILING")
endif()

list(APPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/deps/bison/bison-install")
list(APPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/deps/flex/flex-install")
if (APPLE)
  list(APPEND CMAKE_PREFIX_PATH "/usr/local/opt/bison")
  list(APPEND CMAKE_PREFIX_PATH "/usr/local/opt/flex")
endif()
find_package(BISON 3.4.2 REQUIRED)
find_package(FLEX 2.6.4 REQUIRED)

if (BISON_FOUND)
  get_filename_component(BISON_PARENT_DIR "${BISON_EXECUTABLE}" DIRECTORY)
  message("-- Adding bison lib: ${BISON_PARENT_DIR}/../lib")
  link_directories("${BISON_PARENT_DIR}/../lib/")
endif()

find_library(LIBRT rt)

message("-- FOUND FLEX EXECUTABLE: ${FLEX_EXECUTABLE}")
message("-- FOUND FLEX INCLUDE DIRS: ${FLEX_INCLUDE_DIRS}")

# TODO: Use find_package
# find package
#list(APPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/smt-switch")

set(SMT_SWITCH_DIR "${PROJECT_SOURCE_DIR}/deps/smt-switch")

# rely on cmake modules from smt-switch (required anyway)
set(CMAKE_MODULE_PATH ${SMT_SWITCH_DIR}/cmake)

find_package(GMP REQUIRED)

# Check that dependencies are there
if (NOT EXISTS "${SMT_SWITCH_DIR}/local/include/smt-switch/smt.h")
  message(FATAL_ERROR "Missing smt-switch headers -- try running ./contrib/setup-smt-switch.sh")
endif()

if (NOT EXISTS "${SMT_SWITCH_DIR}/local/lib/libsmt-switch.a")
  message(FATAL_ERROR "Missing smt-switch library -- try running ./contrib/setup-smt-switch.sh")
endif()

if (NOT EXISTS "${SMT_SWITCH_DIR}/local/lib/libsmt-switch-btor.a")
  message(FATAL_ERROR "Missing smt-switch boolector library -- try running ./contrib/setup-smt-switch.sh")
endif()

if (NOT EXISTS "${SMT_SWITCH_DIR}/local/lib/libsmt-switch-cvc5.a")
  message(FATAL_ERROR "Missing smt-switch cvc5 library -- try running ./contrib/setup-smt-switch.sh --with-cvc5")
endif()

if (WITH_MSAT)
  if (NOT EXISTS "${SMT_SWITCH_DIR}/local/lib/libsmt-switch-msat.a")
    message(FATAL_ERROR "Missing smt-switch mathsat library -- try running ./contrib/setup-smt-switch.sh --with-msat")
  endif()
endif()

if (WITH_MSAT_IC3IA)
  if (NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/ic3ia/build/libic3ia.a")
    message(FATAL_ERROR "Missing ic3ia library -- try running ./contrib/setup-ic3ia.sh")
  endif()
endif()

if (WITH_COREIR)
  if (NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/coreir/local/lib/libcoreir.${SHARED_LIB_EXT}")
    message(FATAL_ERROR "Missing coreir library. Try running ./contrib/setup-coreir.sh")
  endif()
  if (NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/coreir/local/lib/libverilogAST.${SHARED_LIB_EXT}")
    message(FATAL_ERROR "Missing verilogAST library from coreir. Try running ./contrib/setup-coreir.sh")
  endif()
endif()


if (NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/btor2tools/src")
  message(FATAL_ERROR "Missing btor2tools -- try running ./contrib/setup-btor2tools.sh")
endif()

# Generate the SMV parser using Flex/Bison
BISON_TARGET(SMVParser ${PROJECT_SOURCE_DIR}/frontends/smvparser.y
  ${CMAKE_CURRENT_BINARY_DIR}/smvparser.cpp
  DEFINES_FILE ${CMAKE_CURRENT_BINARY_DIR}/smvparser.h)

FLEX_TARGET(SMVScanner ${PROJECT_SOURCE_DIR}/frontends/smvparser.l
  ${CMAKE_CURRENT_BINARY_DIR}/smvlexer.cpp)

ADD_FLEX_BISON_DEPENDENCY(SMVScanner SMVParser)

set(INCLUDE_DIRS
  "${PROJECT_SOURCE_DIR}"
  "${PROJECT_SOURCE_DIR}/contrib/"
  "${PROJECT_SOURCE_DIR}/contrib/optionparser-1.7/src"
  "${PROJECT_SOURCE_DIR}/deps/btor2tools/src"
  "${SMT_SWITCH_DIR}/local/include"
  # generated Bison headers go in build directory
  "${CMAKE_CURRENT_BINARY_DIR}"
  "${GMP_INCLUDE_DIR}"
  "${FLEX_INCLUDE_DIRS}"
  )

set(SOURCES
  "${PROJECT_SOURCE_DIR}/core/ts.cpp"
  "${PROJECT_SOURCE_DIR}/core/rts.cpp"
  "${PROJECT_SOURCE_DIR}/core/fts.cpp"
  "${PROJECT_SOURCE_DIR}/core/unroller.cpp"
  "${PROJECT_SOURCE_DIR}/core/functional_unroller.cpp"
  "${PROJECT_SOURCE_DIR}/core/proverresult.cpp"
  "${PROJECT_SOURCE_DIR}/engines/prover.cpp"
  "${PROJECT_SOURCE_DIR}/engines/bmc.cpp"
  "${PROJECT_SOURCE_DIR}/engines/bmc_simplepath.cpp"
  "${PROJECT_SOURCE_DIR}/engines/cegar_ops_uf.cpp"
  "${PROJECT_SOURCE_DIR}/engines/cegar_values.cpp"
  "${PROJECT_SOURCE_DIR}/engines/ceg_prophecy_arrays.cpp"
  "${PROJECT_SOURCE_DIR}/engines/ic3.cpp"
  "${PROJECT_SOURCE_DIR}/engines/ic3base.cpp"
  "${PROJECT_SOURCE_DIR}/engines/ic3bits.cpp"
  "${PROJECT_SOURCE_DIR}/engines/ic3ia.cpp"
  "${PROJECT_SOURCE_DIR}/engines/ic3sa.cpp"
  "${PROJECT_SOURCE_DIR}/engines/interpolantmc.cpp"
  "${PROJECT_SOURCE_DIR}/engines/kinduction.cpp"
  "${PROJECT_SOURCE_DIR}/engines/mbic3.cpp"
  "${PROJECT_SOURCE_DIR}/engines/syguspdr.cpp"
  "${PROJECT_SOURCE_DIR}/frontends/btor2_encoder.cpp"
  "${PROJECT_SOURCE_DIR}/frontends/smv_encoder.cpp"
  "${PROJECT_SOURCE_DIR}/frontends/smv_node.cpp"
  "${PROJECT_SOURCE_DIR}/frontends/vmt_encoder.cpp"
  "${PROJECT_SOURCE_DIR}/modifiers/array_abstractor.cpp"
  "${PROJECT_SOURCE_DIR}/modifiers/control_signals.cpp"
  "${PROJECT_SOURCE_DIR}/modifiers/implicit_predicate_abstractor.cpp"
  "${PROJECT_SOURCE_DIR}/modifiers/history_modifier.cpp"
  "${PROJECT_SOURCE_DIR}/modifiers/mod_ts_prop.cpp"
  "${PROJECT_SOURCE_DIR}/modifiers/ops_abstractor.cpp"
  "${PROJECT_SOURCE_DIR}/modifiers/prophecy_modifier.cpp"
  "${PROJECT_SOURCE_DIR}/modifiers/static_coi.cpp"
  "${PROJECT_SOURCE_DIR}/modifiers/op_abstractor.cpp"
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

if (WITH_COREIR OR WITH_COREIR_EXTERN)
  add_definitions(-DWITH_COREIR)
  set(SOURCES "${SOURCES}" "${PROJECT_SOURCE_DIR}/frontends/coreir_encoder.cpp")
endif()

if (WITH_MSAT)
  # need direct access to mathsat for setting options
  set(INCLUDE_DIRS "${INCLUDE_DIRS}" "${PROJECT_SOURCE_DIR}/deps/mathsat/include/")
endif()

if (WITH_MSAT_IC3IA)
  add_definitions(-DWITH_MSAT_IC3IA)
  # can't include ic3ia directly because some names clash e.g. ic3.h
  set(INCLUDE_DIRS "${INCLUDE_DIRS}" "${PROJECT_SOURCE_DIR}/deps/")
  set(INCLUDE_DIRS "${INCLUDE_DIRS}" "${PROJECT_SOURCE_DIR}/deps/mathsat/include/")
  set(SOURCES "${SOURCES}" "${PROJECT_SOURCE_DIR}/engines/msat_ic3ia.cpp")
endif()


add_library(pono-lib "${PONO_LIB_TYPE}" ${SOURCES})
set_target_properties(pono-lib PROPERTIES OUTPUT_NAME pono)

target_include_directories(pono-lib PUBLIC ${INCLUDE_DIRS})

# Should we build python bindings
option (BUILD_PYTHON_BINDINGS
   "Build Python bindings")

if (BUILD_PYTHON_BINDINGS)
  add_subdirectory(python)
endif()

target_link_libraries(pono-lib PUBLIC "${SMT_SWITCH_DIR}/local/lib/libsmt-switch-btor.a")
target_link_libraries(pono-lib PUBLIC "${SMT_SWITCH_DIR}/local/lib/libsmt-switch-cvc5.a")

if (WITH_MSAT)
  target_link_libraries(pono-lib PUBLIC "${SMT_SWITCH_DIR}/local/lib/libsmt-switch-msat.a")
endif()

if (WITH_MSAT_IC3IA)
  target_link_libraries(pono-lib PUBLIC "${PROJECT_SOURCE_DIR}/deps/ic3ia/build/libic3ia.a")
endif()

if (WITH_COREIR_EXTERN)
  if(APPLE)
    set(COREIR_EXTERN_PREFIX "/usr/local")
  else()
    set(COREIR_EXTERN_PREFIX "/usr")
  endif()
  add_library(coreir SHARED IMPORTED)
  set_target_properties(coreir PROPERTIES
    IMPORTED_LOCATION "${COREIR_EXTERN_PREFIX}/lib/libcoreir${CMAKE_SHARED_LIBRARY_SUFFIX}"
    INTERFACE_INCLUDE_DIRECTORIES "${COREIR_EXTERN_PREFIX}/include"
  )
  add_library(verilogAST SHARED IMPORTED)
  set_target_properties(verilogAST PROPERTIES
    IMPORTED_LOCATION "${COREIR_EXTERN_PREFIX}/lib/libverilogAST${CMAKE_SHARED_LIBRARY_SUFFIX}"
    INTERFACE_INCLUDE_DIRECTORIES "${COREIR_EXTERN_PREFIX}/include"
  )
  target_link_libraries(pono-lib PUBLIC coreir verilogAST)
endif()

if (WITH_COREIR)
  target_include_directories(pono-lib PUBLIC "${PROJECT_SOURCE_DIR}/deps/coreir/local/include")
  target_link_libraries(pono-lib PUBLIC "${PROJECT_SOURCE_DIR}/deps/coreir/local/lib/libcoreir.${SHARED_LIB_EXT}")
  target_link_libraries(pono-lib PUBLIC "${PROJECT_SOURCE_DIR}/deps/coreir/local/lib/libverilogAST.${SHARED_LIB_EXT}")
endif()

target_link_libraries(pono-lib PUBLIC "${SMT_SWITCH_DIR}/local/lib/libsmt-switch.a")
target_link_libraries(pono-lib PUBLIC "${PROJECT_SOURCE_DIR}/deps/btor2tools/build/lib/libbtor2parser.a")
target_link_libraries(pono-lib PUBLIC "${GMPXX_LIBRARIES}")
target_link_libraries(pono-lib PUBLIC "${GMP_LIBRARIES}")
target_link_libraries(pono-lib PUBLIC pthread)
target_link_libraries(pono-lib PUBLIC y)

if (LIBRT)
  target_link_libraries(pono-lib PUBLIC ${LIBRT})
endif()

if (GOOGLE_PERF)
  target_link_libraries(pono-lib PUBLIC ${GOOGLE_PERF})
endif()

enable_testing()
# Add tests subdirectory
# The CMakeLists.txt file there sets up googletest
# and builds all the parametrized tests
add_subdirectory(tests)

add_executable(pono-bin "${PROJECT_SOURCE_DIR}/pono.cpp")
set_target_properties(pono-bin PROPERTIES OUTPUT_NAME pono)

target_include_directories(pono-bin PUBLIC
  "${PROJECT_SOURCE_DIR}/utils"
  "${PROJECT_SOURCE_DIR}/core"
  "${PROJECT_SOURCE_DIR}/contrib/"
  "${PROJECT_SOURCE_DIR}/contrib/optionparser-1.7/src"
  "${PROJECT_SOURCE_DIR}/deps/btor2tools/src"
  "${SMT_SWITCH_DIR}/local/include")

target_link_libraries(pono-bin PUBLIC pono-lib)

if ("${PONO_STATIC_EXEC}" STREQUAL "YES")
  target_link_libraries(pono-bin PUBLIC -static)
endif()

# install smt-switch
install(TARGETS pono-lib DESTINATION lib)
install(TARGETS pono-bin DESTINATION bin)

# install public headers
install(DIRECTORY "${PROJECT_SOURCE_DIR}/core/"
        DESTINATION include/pono/core
        FILES_MATCHING PATTERN "*.h")

install(DIRECTORY "${PROJECT_SOURCE_DIR}/engines/"
        DESTINATION include/pono/engines
        FILES_MATCHING PATTERN "*.h")

install(DIRECTORY "${PROJECT_SOURCE_DIR}/frontends/"
        DESTINATION include/pono/frontends
        FILES_MATCHING PATTERN "*.h")

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/location.hh"
        DESTINATION include/pono/)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/modifiers/"
        DESTINATION include/pono/modifiers
        FILES_MATCHING PATTERN "*.h")

install(DIRECTORY "${PROJECT_SOURCE_DIR}/options/"
        DESTINATION include/pono/options
        FILES_MATCHING PATTERN "*.h")

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/smvparser.h"
        DESTINATION include/pono/)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/utils/"
        DESTINATION include/pono/utils
        FILES_MATCHING PATTERN "*.h")

install(DIRECTORY "${PROJECT_SOURCE_DIR}/contrib/fmt/"
        DESTINATION include/pono/fmt
        FILES_MATCHING PATTERN "*.h")

# uninstall target
# copied from https://gitlab.kitware.com/cmake/community/wikis/FAQ#can-i-do-make-uninstall-with-cmake
if(NOT TARGET uninstall)
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)

    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()
