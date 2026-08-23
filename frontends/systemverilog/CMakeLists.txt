include(FetchContent)

# GIT_SHALLOW is intentionally omitted: FetchContent's shallow clone only
# fetches the tip of each branch, so it cannot reach GIT_TAG once that
# commit is no longer a branch tip (as happens once upstream moves on).
FetchContent_Declare(
  slang
  GIT_REPOSITORY https://github.com/MikePopoloski/slang.git
  GIT_TAG 2e7495bcb64180553ae84d1d78307090823c794c
)
set(SLANG_INCLUDE_TESTS OFF CACHE BOOL "" FORCE)
set(SLANG_INCLUDE_TOOLS OFF CACHE BOOL "" FORCE)
set(SLANG_USE_MIMALLOC OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(slang)
message(STATUS "slang fetched via FetchContent")

add_library(
  pono-systemverilog
  STATIC
  encoder.cpp
  ast_helpers.cpp
  prescan.cpp
  declare.cpp
  instance.cpp
  statement.cpp
  expr.cpp
  sva.cpp
  terms.cpp
)
target_include_directories(pono-systemverilog PUBLIC ${INCLUDE_DIRS})
# slang requires C++20 features in its headers; keep this PRIVATE so it
# doesn't force pono-lib (or pono-bin, or tests) to build at C++20 --
# nothing outside this target includes slang headers directly. (The
# CXX_STANDARD target property is not honored by the Ninja/Makefiles
# generators, so the standard flag is set directly via COMPILE_OPTIONS.)
target_compile_options(pono-systemverilog PRIVATE "${CMAKE_CXX20_STANDARD_COMPILE_OPTION}")
target_link_libraries(pono-systemverilog PUBLIC slang::slang)
