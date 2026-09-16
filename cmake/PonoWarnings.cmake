# Warning policy for pono's own targets, carried by an INTERFACE target that
# each of them links privately.

include(CheckCXXCompilerFlag)

add_library(pono-warnings INTERFACE)

# Adds the flag to pono-warnings if this compiler accepts it, which lets one
# list serve both GCC and Clang: a spelling the compiler does not know is
# skipped. Only positive forms can be probed this way, because both compilers
# accept an unknown -Wno-* silently.
function(pono_add_warning_flag flag)
  # Probe with -Wall on, as the real build has it: some warnings are inert
  # without the group that implies them, and GCC reports those as ignored,
  # which the probe would read as unsupported (-Wformat-security is one).
  set(CMAKE_REQUIRED_FLAGS "-Wall")
  # Clang merely warns about an unknown -W flag, so without this every probe
  # would report success. GCC rejects unknown -W flags on its own and does not
  # recognize this option.
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    string(APPEND CMAKE_REQUIRED_FLAGS " -Werror=unknown-warning-option")
  endif()
  string(MAKE_C_IDENTIFIER "pono_accepts_${flag}" have_flag)
  check_cxx_compiler_flag("${flag}" "${have_flag}")
  if(${have_flag})
    target_compile_options(pono-warnings INTERFACE "${flag}")
  endif()
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

foreach(flag IN LISTS pono_warnings)
  pono_add_warning_flag("${flag}")
endforeach()

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

foreach(name IN LISTS pono_fatal_warnings)
  pono_add_warning_flag("-Werror=${name}")
endforeach()

# Off by default, so that warnings a future compiler adds cannot break an
# end-user build.
option(PONO_STRICT "Treat compiler warnings as errors")

if(PONO_STRICT)
  pono_add_warning_flag(-Werror)
endif()
