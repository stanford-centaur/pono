cmake_minimum_required(VERSION 3.1)

project(cosa2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/core")
include_directories("${PROJECT_SOURCE_DIR}/frontends")

# TODO: Use find_package
# find package
#list(APPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/smt-switch")

add_executable(cosa2
  "${PROJECT_SOURCE_DIR}/cosa2.cpp"
  "${PROJECT_SOURCE_DIR}/btor2tools/build/libbtor2parser.a"
  "${PROJECT_SOURCE_DIR}/frontends/btor2_encoder.cpp"
  "${PROJECT_SOURCE_DIR}/core/rts.cpp"
  "${PROJECT_SOURCE_DIR}/core/fts.cpp"
  "${PROJECT_SOURCE_DIR}/unroller.cpp"
  "${PROJECT_SOURCE_DIR}/bmc.cpp"
  "${PROJECT_SOURCE_DIR}/kinduction.cpp"
  "${PROJECT_SOURCE_DIR}/core/prop.cpp"
  )

target_include_directories(cosa2 PUBLIC "${PROJECT_SOURCE_DIR}/utils")
target_include_directories(cosa2 PUBLIC "${PROJECT_SOURCE_DIR}/core")
target_include_directories(cosa2 PUBLIC "${PROJECT_SOURCE_DIR}/btor2tools/src")
target_include_directories(cosa2 PUBLIC "${PROJECT_SOURCE_DIR}/smt-switch/local/include")
target_include_directories(cosa2 PUBLIC "${PROJECT_SOURCE_DIR}/contrib/optionparser-1.7/src")

target_link_directories(cosa2 PUBLIC "${PROJECT_SOURCE_DIR}/smt-switch/local/lib")
target_link_libraries(cosa2 smt-switch)
target_link_libraries(cosa2 smt-switch-btor)
