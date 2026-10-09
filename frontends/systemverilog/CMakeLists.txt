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
# tomlplusplus is header-only by default, but a distro may build it in
# its compiled mode instead -- Arch ships libtomlplusplus.so and no
# archive -- and slang's FIND_PACKAGE_ARGS would then take a shared
# object a fully static pono cannot link. Nothing recommends the system
# copy here the way an archived fmt is worth keeping: the fetched
# tomlplusplus is an INTERFACE target, so it costs a clone and no link
# at all.
if(PONO_LINK_STATIC)
  set(CMAKE_DISABLE_FIND_PACKAGE_tomlplusplus ON)
endif()
FetchContent_MakeAvailable(slang)
message(STATUS "slang fetched via FetchContent")

# An OBJECT library rather than a STATIC one: these objects belong to
# pono-lib and are linked straight into it. As a separate archive they
# were mutually dependent with it -- calling into TransitionSystem
# while pono-lib linked them back in -- which a linker taking one pass
# over each archive can only resolve if it happens to reach them in
# the right order. CMake will not let the two declare that cycle,
# since it is allowed only between static libraries and pono-lib is
# shared by default.
add_library(
  pono-systemverilog
  OBJECT
  encoder.cpp
  ast_helpers.cpp
  bit_utils.cpp
  symbol_table.cpp
  declarer.cpp
  instance_encoder.cpp
  statement_encoder.cpp
  expr_encoder.cpp
  assertion_walker.cpp
  tableau.cpp
)
target_include_directories(pono-systemverilog PUBLIC ${INCLUDE_DIRS})
# slang requires C++20 features in its headers; keep this PRIVATE so it
# doesn't force pono-lib (or pono-bin, or tests) to build at C++20 --
# nothing outside this target includes slang headers directly. (The
# CXX_STANDARD target property is not honored by the Ninja/Makefiles
# generators, so the standard flag is set directly via COMPILE_OPTIONS.)
target_compile_options(pono-systemverilog PRIVATE "${CMAKE_CXX20_STANDARD_COMPILE_OPTION}")
target_link_libraries(pono-systemverilog PUBLIC slang::slang smt-switch fmt::fmt-header-only)
