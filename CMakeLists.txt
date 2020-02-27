cmake_minimum_required(VERSION 3.1)

project(cosa2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (WITH_MSAT)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_MSAT")
endif()

include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/core")
include_directories("${PROJECT_SOURCE_DIR}/engines")
include_directories("${PROJECT_SOURCE_DIR}/frontends")
include_directories("${PROJECT_SOURCE_DIR}/options")
include_directories("/usr/local/include")

link_directories("/usr/local/lib")

# TODO: Use find_package
# find package
#list(APPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/smt-switch")

add_library(cosa2
  "${PROJECT_SOURCE_DIR}/frontends/btor2_encoder.cpp"
  "${PROJECT_SOURCE_DIR}/core/ts.cpp"
  "${PROJECT_SOURCE_DIR}/core/rts.cpp"
  "${PROJECT_SOURCE_DIR}/core/fts.cpp"
  "${PROJECT_SOURCE_DIR}/core/prop.cpp"
  "${PROJECT_SOURCE_DIR}/core/unroller.cpp"
  "${PROJECT_SOURCE_DIR}/engines/prover.cpp"
  "${PROJECT_SOURCE_DIR}/engines/bmc.cpp"
  "${PROJECT_SOURCE_DIR}/engines/interpolant.cpp"
  "${PROJECT_SOURCE_DIR}/engines/bmc_simplepath.cpp"
  "${PROJECT_SOURCE_DIR}/engines/kinduction.cpp"
  "${PROJECT_SOURCE_DIR}/utils/logger.cpp"
  "${PROJECT_SOURCE_DIR}/utils/term_analysis.cpp"
  )

if (APPLE)
  set(SHARED_LIBRARY_EXT dylib)
else()
  set(SHARED_LIBRARY_EXT so)
endif()

# Check that dependencies are there
if (NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/include/smt-switch/smt.h")
  message(FATAL_ERROR "Missing smt-switch headers -- try running ./contrib/setup-smt-switch.sh")
endif()

if (NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch.${SHARED_LIBRARY_EXT}")
  message(FATAL_ERROR "Missing smt-switch library -- try running ./contrib/setup-smt-switch.sh")
endif()

if (NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch-btor.${SHARED_LIBRARY_EXT}")
  message(FATAL_ERROR "Missing smt-switch boolector library -- try running ./contrib/setup-smt-switch.sh")
endif()

if (WITH_MSAT)
  if (NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch-msat.${SHARED_LIBRARY_EXT}")
    message(FATAL_ERROR "Missing smt-switch mathsat library -- try running ./contrib/setup-smt-switch.sh --with-msat")
  endif()
endif()

if (NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/btor2tools/src")
  message(FATAL_ERROR "Missing btor2tools -- try running ./contrib/setup-btor2tools.sh")
endif()

target_include_directories(cosa2 PUBLIC
  "${PROJECT_SOURCE_DIR}/utils"
  "${PROJECT_SOURCE_DIR}/core"
  "${PROJECT_SOURCE_DIR}/contrib/"
  "${PROJECT_SOURCE_DIR}/contrib/optionparser-1.7/src"
  "${PROJECT_SOURCE_DIR}/deps/btor2tools/src"
  "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/include")

target_link_libraries(cosa2 PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch.${SHARED_LIBRARY_EXT}")
target_link_libraries(cosa2 PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch-btor.${SHARED_LIBRARY_EXT}")
if (WITH_MSAT)
  target_link_libraries(cosa2 PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch-msat.${SHARED_LIBRARY_EXT}")
endif()
target_link_libraries(cosa2 PUBLIC "${PROJECT_SOURCE_DIR}/deps/btor2tools/build/lib/libbtor2parser.a")
target_link_libraries(cosa2 PUBLIC gmp)

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

target_link_libraries(cosa2-bin PUBLIC cosa2)
target_link_libraries(cosa2-bin PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch.${SHARED_LIBRARY_EXT}")
target_link_libraries(cosa2-bin PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch-btor.${SHARED_LIBRARY_EXT}")
if (WITH_MSAT)
  target_link_libraries(cosa2-bin PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch-msat.${SHARED_LIBRARY_EXT}")
endif()
target_link_libraries(cosa2-bin PUBLIC "${PROJECT_SOURCE_DIR}/deps/btor2tools/build/lib/libbtor2parser.a")
target_link_libraries(cosa2-bin PUBLIC gmp)
