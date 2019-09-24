cmake_minimum_required(VERSION 3.1)

project(cosa2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/core")
include_directories("${PROJECT_SOURCE_DIR}/engines")
include_directories("${PROJECT_SOURCE_DIR}/frontends")
include_directories("${PROJECT_SOURCE_DIR}/options")

# TODO: Use find_package
# find package
#list(APPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/smt-switch")

add_executable(cosa2
  "${PROJECT_SOURCE_DIR}/cosa2.cpp"
  "${PROJECT_SOURCE_DIR}/deps/btor2tools/build/libbtor2parser.a"
  "${PROJECT_SOURCE_DIR}/frontends/btor2_encoder.cpp"
  "${PROJECT_SOURCE_DIR}/core/rts.cpp"
  "${PROJECT_SOURCE_DIR}/core/fts.cpp"
  "${PROJECT_SOURCE_DIR}/core/prop.cpp"
  "${PROJECT_SOURCE_DIR}/core/unroller.cpp"
  "${PROJECT_SOURCE_DIR}/engines/bmc.cpp"
  "${PROJECT_SOURCE_DIR}/engines/kinduction.cpp"
  )

# Check that dependencies are there
if (NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/include/smt-switch/smt.h")
  message(FATAL_ERROR "Missing smt-switch headers -- try running ./contrib/setup-smt-switch.sh")
endif()

if (NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch.so")
  message(FATAL_ERROR "Missing smt-switch library -- try running ./contrib/setup-smt-switch.sh")
endif()

if (NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch-btor.so")
  message(FATAL_ERROR "Missing smt-switch boolector library -- try running ./contrib/setup-smt-switch.sh")
endif()

if (NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/btor2tools/src")
  message(FATAL_ERROR "Missing btor2tools -- try running ./contrib/setup-btor2tools.sh")
endif()

target_include_directories(cosa2 PUBLIC "${PROJECT_SOURCE_DIR}/utils")
target_include_directories(cosa2 PUBLIC "${PROJECT_SOURCE_DIR}/core")
target_include_directories(cosa2 PUBLIC "${PROJECT_SOURCE_DIR}/contrib/optionparser-1.7/src")
target_include_directories(cosa2 PUBLIC "${PROJECT_SOURCE_DIR}/deps/btor2tools/src")
target_include_directories(cosa2 PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/include")

target_link_libraries(cosa2 PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch.so")
target_link_libraries(cosa2 PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/lib/libsmt-switch-btor.so")
