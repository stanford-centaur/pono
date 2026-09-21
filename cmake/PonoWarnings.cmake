# Warning policy for pono's own targets, carried by an INTERFACE target that
# each of them links privately.

add_library(pono-warnings INTERFACE)

# Off by default, so that warnings a future compiler adds cannot break an
# end-user build.
option(PONO_STRICT "Treat compiler warnings as errors")

# Sets <out_var> to the flags among its arguments that this compiler accepts,
# which lets one list serve both GCC and Clang: a spelling the compiler does
# not know is skipped. Only positive forms can be filtered this way, because
# both compilers accept an unknown -Wno-* silently.
#
# One compiler run per flag answers this. Recognizing an option needs no link
# step and no headers, so a syntax-only pass over a trivial source is enough.
function(pono_filter_supported_flags out_var)
  # CXX="ccache g++" reaches us split in two, and the first half alone is not
  # a compiler.
  separate_arguments(argv1 NATIVE_COMMAND "${CMAKE_CXX_COMPILER_ARG1}")
  set(compiler "${CMAKE_CXX_COMPILER}" ${argv1})

  # Probe with -Wall on, as the real build has it: some warnings are inert
  # without the group that implies them, and GCC reports those as ignored,
  # which the probe would read as unsupported (-Wformat-security is one).
  set(context -Wall)
  # Clang merely warns about an unknown -W flag, so without this every probe
  # would report success. GCC rejects unknown -W flags on its own and does not
  # recognize this option.
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    list(APPEND context -Werror=unknown-warning-option)
  endif()

  set(source "${CMAKE_BINARY_DIR}/CMakeFiles/pono-flag-probe.cpp")
  file(WRITE "${source}" "int main(void) { return 0; }\n")

  # A compiler that cannot run here at all rejects every flag, which looks
  # exactly like a compiler that supports none of them. Say so, rather than
  # building unwarned in silence.
  execute_process(
    COMMAND ${compiler} ${context} -fsyntax-only "${source}"
    RESULT_VARIABLE status
    OUTPUT_QUIET
    ERROR_QUIET
  )
  if(NOT status EQUAL 0)
    message(
      WARNING
      "Cannot run ${CMAKE_CXX_COMPILER} to find out which warning flags it "
      "takes, so pono is built without its warning set."
    )
    set(${out_var} "" PARENT_SCOPE)
    return()
  endif()

  set(supported "")
  foreach(flag IN LISTS ARGN)
    execute_process(
      COMMAND ${compiler} ${context} "${flag}" -fsyntax-only "${source}"
      RESULT_VARIABLE status
      OUTPUT_QUIET
      ERROR_QUIET
    )
    if(status EQUAL 0)
      list(APPEND supported "${flag}")
    endif()
  endforeach()

  set(${out_var} "${supported}" PARENT_SCOPE)
endfunction()

# The broad groups, then individual warnings outside them that have caught real
# bugs here or that guard against the classes those bugs came from.
set(
  pono_warnings
  -Wall # the usual baseline
  -Wextra # and the tier above it
  -Wpedantic # constructs outside strict ISO C++
  -Wc++20-compat # code whose meaning would change in C++20
  -Wcast-align # a cast that raises the required alignment
  -Wcast-qual # a cast that drops const or volatile
  -Wconversion # an implicit conversion that can lose value
  -Wdelete-non-abstract-non-virtual-dtor # non-final, non-virtual dtor
  -Wdeprecated-copy-dtor # implicit copy where a dtor is declared
  -Wdeprecated-copy-with-dtor # the same, as Clang spells it
  -Wdeprecated-dynamic-exception-spec # throw() rather than noexcept
  -Wduplicated-branches # if and else with identical bodies
  -Wduplicated-cond # the same condition twice in one if chain
  -Wextra-semi # a stray ; after a member function body
  -Wextra-semi-stmt # a ; that forms an empty statement
  -Wformat=2 # printf arguments, including non-literal formats
  -Wheader-hygiene # using namespace at header scope
  -Wimplicit-int-float-conversion # an int to float losing precision
  -Wlogical-op # && or || where & or | was likely meant
  -Wnon-virtual-dtor # a polymorphic class without a virtual dtor
  -Wnull-dereference # a path that dereferences null
  -Woverloaded-virtual # hides a base overload instead of overriding
  -Wshadow # a declaration that hides an outer one
  -Wshadow-field # a parameter or local that hides a member
  -Wshorten-64-to-32 # a 64-bit value narrowed to 32
  -Wsign-conversion # signedness changed by an implicit conversion
  -Wsuggest-destructor-override # a dtor overriding without saying so
  -Wsuggest-override # the same for any member function
  -Wundef # #if on a macro that is not defined
  -Wunreachable-code-break # a break that cannot be reached
  -Wunused-exception-parameter # a caught exception never used
  -Wvla # a variable-length array
)

# Fatal in every build, not only under PONO_STRICT. Each of these marks a
# construct that is a bug wherever it appears, whose consequence is silent
# misbehavior or undefined behavior rather than style.
set(
  pono_fatal_warnings
  dangling-gsl # a pointer kept to a temporary that has died
  dangling-pointer # the same, as GCC spells it
  delete-non-virtual-dtor # delete through a base with a non-virtual dtor
  div-by-zero # constant division by zero
  division-by-zero # the same, as Clang spells it
  dynamic-class-memaccess # memcpy over an object, clobbering its vptr
  format-security # printf given a format string it cannot check
  free-nonheap-object # free() of something never allocated
  inconsistent-missing-override # overrides a virtual without saying so
  infinite-recursion # a function that can only call itself
  int-in-bool-context # arithmetic where a condition was meant
  memset-transposed-args # memset(p, size, value), the last two swapped
  mismatched-new-delete # new[] released with a plain delete
  nonnull # null passed where the callee declares it cannot be
  placement-new # placement new into too small a buffer
  return-local-addr # returns the address of a local
  return-stack-address # the same, as Clang spells it
  return-type # control reaches the end of a non-void function
  self-assign # x = x
  self-move # x = std::move(x)
  sequence-point # i = i++, evaluation order undefined
  sizeof-array-argument # sizeof a parameter that decayed to a pointer
  sizeof-pointer-memaccess # memcpy(p, q, sizeof(p)), the pointer's size
  string-compare # compares literal addresses rather than contents
  switch # an enum switch missing a case
  tautological-compare # a comparison whose result cannot vary
  uninitialized # read before anything is written
  unsequenced # the same as sequence-point, as Clang spells it
  unused-result # discards a value the callee marked must-use
  user-defined-warnings # a diagnostic the code itself raises
  varargs # va_arg used with a type it cannot receive
)

list(TRANSFORM pono_fatal_warnings PREPEND "-Werror=" OUTPUT_VARIABLE pono_fatal_flags)
set(pono_candidate_flags ${pono_warnings} ${pono_fatal_flags})

# Filtering spends a compiler run per flag, so the answer is cached and only
# redone when something it depends on changes. The compiler's size and
# timestamp are part of that: CMake keeps the version it detected when the
# cache was first written, so an in-place upgrade would otherwise go
# unnoticed. Clear PONO_WARNING_FLAGS_SIGNATURE to force a fresh look.
get_filename_component(pono_cxx_path "${CMAKE_CXX_COMPILER}" REALPATH)
file(TIMESTAMP "${pono_cxx_path}" pono_cxx_mtime UTC)
file(SIZE "${pono_cxx_path}" pono_cxx_size)
string(
  SHA256 pono_warning_signature
  "${CMAKE_CXX_COMPILER}|${CMAKE_CXX_COMPILER_ARG1}|${CMAKE_CXX_COMPILER_ID}\
|${CMAKE_CXX_COMPILER_VERSION}|${pono_cxx_mtime}|${pono_cxx_size}\
|${pono_candidate_flags}"
)

# Accepting none of them is a legitimate answer, so what decides whether the
# cache can be reused is that both entries exist, not that they are non-empty.
if(
  NOT DEFINED PONO_WARNING_FLAGS
  OR NOT DEFINED PONO_WARNING_FLAGS_SIGNATURE
  OR NOT PONO_WARNING_FLAGS_SIGNATURE STREQUAL pono_warning_signature
)
  # Every flag in the lists above is a GCC or Clang spelling, as is the
  # -fsyntax-only the filter probes with.
  if(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|.*Clang)$")
    pono_filter_supported_flags(pono_supported_flags ${pono_candidate_flags})
  else()
    set(pono_supported_flags "")
  endif()
  set(
    PONO_WARNING_FLAGS
    "${pono_supported_flags}"
    CACHE INTERNAL
    "Warning flags this compiler accepts, of the ones pono asks for"
  )
  set(
    PONO_WARNING_FLAGS_SIGNATURE
    "${pono_warning_signature}"
    CACHE INTERNAL
    "Compiler and flag list that PONO_WARNING_FLAGS was determined for"
  )
endif()

target_compile_options(pono-warnings INTERFACE ${PONO_WARNING_FLAGS})

# Set directly, since both compilers accept any -Wno-* silently and so cannot
# be filtered for it, and after the flags above so it wins over the -Wextra
# that turns it on. Every parameter it reports is held in place by a virtual,
# template or callback signature, leaving nothing to act on.
target_compile_options(pono-warnings INTERFACE -Wno-unused-parameter)

# Not filtered: -Werror predates every GCC and Clang that can build pono, and
# leaving PONO_STRICT out of the signature keeps toggling it free.
if(PONO_STRICT)
  target_compile_options(pono-warnings INTERFACE -Werror)
endif()
