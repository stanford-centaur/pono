cmake_minimum_required(VERSION 3.5.1)

project(cosa2)

set(COSA2_MAJOR 0)   # Major component of Cosa2 version
set(COSA2_MINOR 1)   # Minor component of Cosa2 version
set(COSA2_RELEASE 0) # Release component of Cosa2 version

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

add_definitions(-DCOSA2_SRC_DIR=${PROJECT_SOURCE_DIR})

if (WITH_MSAT)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_MSAT")
endif()

include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/utils")
include_directories("${PROJECT_SOURCE_DIR}/core")
include_directories("${PROJECT_SOURCE_DIR}/engines")
include_directories("${PROJECT_SOURCE_DIR}/frontends")
include_directories("${PROJECT_SOURCE_DIR}/options")
include_directories("/usr/local/include")
include_directories("/usr/local/Cellar/gmp/6.2.0/include")

link_directories("/usr/local/lib")
link_directories("/usr/local/Cellar/gmp/6.2.0/lib")

list(APPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/deps/bison/bison-install")
find_package(BISON 3.4.2 REQUIRED)
find_package(FLEX 2.6.4 REQUIRED)

if (BISON_FOUND)
  link_directories("${PROJECT_SOURCE_DIR}/deps/bison/bison-install/lib/")
endif()

# TODO: Use find_package
# find package
#list(APPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/smt-switch")

# Check that dependencies are there
if (NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/include/smt-switch/smt.h")
  message(FATAL_ERROR "Missing smt-switch headers -- try running ./contrib/setup-smt-switch.sh")
endif()

if (NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch.a")
  message(FATAL_ERROR "Missing smt-switch library -- try running ./contrib/setup-smt-switch.sh")
endif()

if (NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch-btor.a")
  message(FATAL_ERROR "Missing smt-switch boolector library -- try running ./contrib/setup-smt-switch.sh")
endif()

if (WITH_MSAT)
  if (NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch-msat.a")
    message(FATAL_ERROR "Missing smt-switch mathsat library -- try running ./contrib/setup-smt-switch.sh --with-msat")
  endif()
endif()

if (WITH_CVC4)
  if (NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch-cvc4.a")
    message(FATAL_ERROR "Missing smt-switch cvc4 library -- try running ./contrib/setup-smt-switch.sh --with-cvc4")
  endif()
endif()

if (WITH_COREIR)
  if (NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/coreir/local/lib/libcoreir.so")
    message(FATAL_ERROR "Missing coreir library. Try running ./contrib/setup-coreir.sh")
  endif()
  if (NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/coreir/local/lib/libverilogAST.so")
    message(FATAL_ERROR "Missing verilogAST library from coreir. Try running ./contrib/setup-coreir.sh")
  endif()
endif()


if (NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/btor2tools/src")
  message(FATAL_ERROR "Missing btor2tools -- try running ./contrib/setup-btor2tools.sh")
endif()

# rely on cmake modules from smt-switch (required anyway)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/deps/smt-switch/cmake)

# Generate the SMV parser using Flex/Bison
BISON_TARGET(SMVParser ${PROJECT_SOURCE_DIR}/frontends/smvparser.y
  ${CMAKE_CURRENT_BINARY_DIR}/smvparser.cpp
  DEFINES_FILE ${CMAKE_CURRENT_BINARY_DIR}/smvparser.h)

FLEX_TARGET(SMVScanner ${PROJECT_SOURCE_DIR}/frontends/smvparser.l
  ${CMAKE_CURRENT_BINARY_DIR}/smvlexer.cpp)

ADD_FLEX_BISON_DEPENDENCY(SMVScanner SMVParser)

add_library(cosa2-lib "${COSA2_LIB_TYPE}"
  "${PROJECT_SOURCE_DIR}/core/ts.cpp"
  "${PROJECT_SOURCE_DIR}/core/rts.cpp"
  "${PROJECT_SOURCE_DIR}/core/fts.cpp"
  "${PROJECT_SOURCE_DIR}/core/prop.cpp"
  "${PROJECT_SOURCE_DIR}/core/unroller.cpp"
  "${PROJECT_SOURCE_DIR}/core/proverresult.cpp"
  "${PROJECT_SOURCE_DIR}/engines/prover.cpp"
  "${PROJECT_SOURCE_DIR}/engines/bmc.cpp"
  "${PROJECT_SOURCE_DIR}/engines/interpolantmc.cpp"
  "${PROJECT_SOURCE_DIR}/engines/bmc_simplepath.cpp"
  "${PROJECT_SOURCE_DIR}/engines/kinduction.cpp"
  "${PROJECT_SOURCE_DIR}/utils/logger.cpp"
  "${PROJECT_SOURCE_DIR}/utils/term_analysis.cpp"
  "${PROJECT_SOURCE_DIR}/printers/vcd_witness_printer.cpp"
  "${PROJECT_SOURCE_DIR}/frontends/btor2_encoder.cpp"
  "${PROJECT_SOURCE_DIR}/frontends/smv_encoder.cpp"
  "${BISON_SMVParser_OUTPUTS}"
  "${FLEX_SMVScanner_OUTPUTS}"
  )

if (WITH_COREIR)
  add_definitions(-DWITH_COREIR)
  target_sources(cosa2-lib PUBLIC "${PROJECT_SOURCE_DIR}/frontends/coreir_encoder.cpp")
endif()

target_include_directories(cosa2-lib PUBLIC
  "${PROJECT_SOURCE_DIR}/utils"
  "${PROJECT_SOURCE_DIR}/core"
  "${PROJECT_SOURCE_DIR}/contrib/"
  "${PROJECT_SOURCE_DIR}/contrib/optionparser-1.7/src"
  "${PROJECT_SOURCE_DIR}/deps/btor2tools/src"
  "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/include"
  # generated Bison headers go in build directory
  "${CMAKE_CURRENT_BINARY_DIR}"
  )

# Should we build python bindings
option (BUILD_PYTHON_BINDINGS
   "Build Python bindings")
option (USE_PYTHON2
   "Use Python2 for bindings")

if (BUILD_PYTHON_BINDINGS)
  if (USE_PYTHON2)
    message("-- We understand: it's hard to let go. We'll try to use Python2, but it's not really supported, so no promises.")
    find_package(PythonInterp 2.7 REQUIRED)
  else()
    find_package(PythonInterp 3 REQUIRED)
  endif()
  add_subdirectory(python)
endif()

target_link_libraries(cosa2-lib PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch-btor.a")

if (WITH_MSAT)
  target_link_libraries(cosa2-lib PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch-msat.a")
endif()

if (WITH_CVC4)
  target_link_libraries(cosa2-lib PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch-cvc4.a")
endif()

if (WITH_COREIR)
  target_include_directories(cosa2-lib PUBLIC "${PROJECT_SOURCE_DIR}/deps/coreir/local/include")
  target_link_libraries(cosa2-lib PUBLIC "${PROJECT_SOURCE_DIR}/deps/coreir/local/lib/libcoreir.so")
  target_link_libraries(cosa2-lib PUBLIC "${PROJECT_SOURCE_DIR}/deps/coreir/local/lib/libverilogAST.so")
endif()

target_link_libraries(cosa2-lib PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch.a")
target_link_libraries(cosa2-lib PUBLIC "${PROJECT_SOURCE_DIR}/deps/btor2tools/build/lib/libbtor2parser.a")
target_link_libraries(cosa2-lib PUBLIC gmp)
target_link_libraries(cosa2-lib PUBLIC pthread)
target_link_libraries(cosa2-lib PUBLIC y)

enable_testing()
# Add tests subdirectory
# The CMakeLists.txt file there sets up googletest
# and builds all the parametrized tests
add_subdirectory(tests)

add_executable(cosa2-bin "${PROJECT_SOURCE_DIR}/cosa2.cpp")
set_target_properties(cosa2-bin PROPERTIES OUTPUT_NAME cosa2)

target_include_directories(cosa2-bin PUBLIC
  "${PROJECT_SOURCE_DIR}/utils"
  "${PROJECT_SOURCE_DIR}/core"
  "${PROJECT_SOURCE_DIR}/contrib/"
  "${PROJECT_SOURCE_DIR}/contrib/optionparser-1.7/src"
  "${PROJECT_SOURCE_DIR}/deps/btor2tools/src"
  "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/include")

target_link_libraries(cosa2-bin PUBLIC cosa2-lib)

if ("${COSA2_STATIC_EXEC}" STREQUAL "YES")
  target_link_libraries(cosa2-bin PUBLIC -static)
endif()

