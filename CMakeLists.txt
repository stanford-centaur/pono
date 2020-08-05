cmake_minimum_required(VERSION 3.1)

project(pono)

set(PONO_MAJOR 0)   # Major component of Pono version
set(PONO_MINOR 1)   # Minor component of Pono version
set(PONO_RELEASE 0) # Release component of Pono version

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

if (WITH_MSAT)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_MSAT")
endif()

include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/utils")
include_directories("${PROJECT_SOURCE_DIR}/core")
include_directories("${PROJECT_SOURCE_DIR}/engines")
include_directories("${PROJECT_SOURCE_DIR}/frontends")
include_directories("${PROJECT_SOURCE_DIR}/smt")
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

message("-- FOUND FLEX EXECUTABLE: ${FLEX_EXECUTABLE}")
message("-- FOUND FLEX INCLUDE DIRS: ${FLEX_INCLUDE_DIRS}")

if (APPLE)
  message("setting include dir")
  include_directories("/usr/local/opt/flex/include")
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

if (NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch-cvc4.a")
  message(FATAL_ERROR "Missing smt-switch cvc4 library -- try running ./contrib/setup-smt-switch.sh --with-cvc4")
endif()

if (WITH_MSAT)
  if (NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch-msat.a")
    message(FATAL_ERROR "Missing smt-switch mathsat library -- try running ./contrib/setup-smt-switch.sh --with-msat")
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

add_library(pono-lib "${PONO_LIB_TYPE}"
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
  "${PROJECT_SOURCE_DIR}/frontends/btor2_encoder.cpp"
  "${PROJECT_SOURCE_DIR}/frontends/smv_encoder.cpp"
  "${PROJECT_SOURCE_DIR}/modifiers/array_abstractor.cpp"
  "${PROJECT_SOURCE_DIR}/printers/vcd_witness_printer.cpp"
  "${PROJECT_SOURCE_DIR}/smt/available_solvers.cpp"
  "${PROJECT_SOURCE_DIR}/utils/logger.cpp"
  "${PROJECT_SOURCE_DIR}/utils/term_analysis.cpp"
  "${PROJECT_SOURCE_DIR}/options/options.cpp"
  "${BISON_SMVParser_OUTPUTS}"
  "${FLEX_SMVScanner_OUTPUTS}"
  )

set_target_properties(pono-lib PROPERTIES OUTPUT_NAME pono)

if (WITH_COREIR)
  add_definitions(-DWITH_COREIR)
  target_sources(pono-lib PUBLIC "${PROJECT_SOURCE_DIR}/frontends/coreir_encoder.cpp")
endif()

target_include_directories(pono-lib PUBLIC
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

target_link_libraries(pono-lib PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch-btor.a")
target_link_libraries(pono-lib PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch-cvc4.a")

if (WITH_MSAT)
  target_link_libraries(pono-lib PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch-msat.a")
endif()

if (WITH_COREIR)
  target_include_directories(pono-lib PUBLIC "${PROJECT_SOURCE_DIR}/deps/coreir/local/include")
  target_link_libraries(pono-lib PUBLIC "${PROJECT_SOURCE_DIR}/deps/coreir/local/lib/libcoreir.so")
  target_link_libraries(pono-lib PUBLIC "${PROJECT_SOURCE_DIR}/deps/coreir/local/lib/libverilogAST.so")
endif()

target_link_libraries(pono-lib PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch.a")
target_link_libraries(pono-lib PUBLIC "${PROJECT_SOURCE_DIR}/deps/btor2tools/build/lib/libbtor2parser.a")
target_link_libraries(pono-lib PUBLIC gmp)
target_link_libraries(pono-lib PUBLIC pthread)
target_link_libraries(pono-lib PUBLIC y)

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
  "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/include")

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

install(DIRECTORY "${PROJECT_SOURCE_DIR}/utils/"
        DESTINATION include/pono/utils
        FILES_MATCHING PATTERN "*.h")

install(DIRECTORY "${PROJECT_SOURCE_DIR}/engines/"
        DESTINATION include/pono/engines
        FILES_MATCHING PATTERN "*.h")

install(DIRECTORY "${PROJECT_SOURCE_DIR}/frontends/"
        DESTINATION include/pono/frontends
        FILES_MATCHING PATTERN "*.h")

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/smvparser.h"
        DESTINATION include/pono/)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/location.hh"
        DESTINATION include/pono/)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/options/"
        DESTINATION include/pono/options
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
