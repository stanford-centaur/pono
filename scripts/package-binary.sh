#!/usr/bin/env bash
# Packages a statically linked pono binary for distribution, for example on
# Zenodo as an SV-COMP tool archive, as a directory pono-<version> and its
# zip archive, where <version> is the git describe output that pono --version
# reports. Beside the binary, the directory holds a LICENSE that appends a
# notice for every third-party component linked into the binary to pono's
# own, the license texts of those components in LICENSES/, and the files in
# scripts/package: README.md filled in from README.md.in, and a smoke test,
# which this script runs on the package.
#
# The components are detected from the symbols of the binary and located
# through the CMake cache of its build. License texts that are not on disk
# are downloaded from upstream at the revision the build used. A linked
# component whose license this script does not cover aborts the packaging.
set -euo pipefail
shopt -s inherit_errexit

usage() {
  cat <<EOF
Usage: $0 [<option> ...] <build-dir>

Packages <build-dir>/pono, which must be statically linked (configure.sh
--static) and not stripped, with the license texts of everything linked
into it.

-h, --help              display this message and exit
-o, --output-dir DIR    write pono-<version>/ and pono-<version>.zip to DIR
                        (default: current directory)
--no-strip              keep the symbol table in the packaged binary
EOF
  exit 0
}

die() {
  echo "*** $0: $*" >&2
  exit 1
}

warn() {
  echo "*** $0: warning: $*" >&2
}

out_dir=.
strip_binary=yes
build_dir=
while [[ $# -gt 0 ]]; do
  case $1 in
    -h | --help) usage ;;
    -o | --output-dir)
      [[ $# -ge 2 ]] || die "missing argument to $1"
      out_dir=$2
      shift
      ;;
    --no-strip) strip_binary=no ;;
    -*) die "unknown option: $1" ;;
    *)
      [[ -z $build_dir ]] || die "unexpected argument: $1"
      build_dir=$1
      ;;
  esac
  shift
done
[[ -n $build_dir ]] || die "missing build directory, see --help"

for tool in curl nm readelf strip zip; do
  command -v "$tool" >/dev/null || die "$tool is required but not installed"
done

# The files that go into the package, as they are or filled in.
assets="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)/package"
cache=$build_dir/CMakeCache.txt
binary=$build_dir/pono
[[ -f $cache ]] || die "$cache not found"
[[ -x $binary ]] || die "$binary not found"

# Prints the value of entry $1 in the CMake cache, or nothing when unset.
cache_entry() {
  sed -n "/^$1:/{s/^[^=]*=//p;q;}" "$cache"
}

# Prints the value of entry $1 in the CMake cache, aborting when unset.
cache_path() {
  local value
  value=$(cache_entry "$1")
  [[ -n $value ]] || die "$1 is not set in $cache"
  echo "$value"
}

# Prints the non-empty value that CMake file $2 sets variable $1 to.
cmake_set() {
  local value
  value=$(sed -n "s/^ *set($1 \"\(..*\)\")\$/\1/p" "$2")
  [[ -n $value ]] || die "$2 sets no $1"
  echo "$value"
}

# Downloads URL $1 to file $2.
download() {
  curl -fsSL --retry 3 -o "$2" "$1" || die "cannot download $1"
}

# Downloads file $3 of GitHub repository $1 at revision $2 to file $4.
github_file() {
  download "https://raw.githubusercontent.com/$1/$2/$3" "$4"
}

program_headers=$(readelf -l "$binary")
if [[ $program_headers == *INTERP* ]]; then
  die "$binary is dynamically linked, configure the build with --static"
fi
build_type=$(cache_entry CMAKE_BUILD_TYPE)
if [[ $build_type != Release ]]; then
  warn "packaging a ${build_type:-default} build, not a Release one"
fi

# pono --version prints "<version> (<git describe output>)", see RELEASING.md,
# and then exits with status 2, as it does after --help.
version_line=$("$binary" --version) || true
version=${version_line%% *}
if [[ ! $version =~ ^[0-9][0-9A-Za-z.+-]*$ ]]; then
  die "cannot read a version from: $version_line"
fi
revision=
if [[ $version_line == *" ("*")" ]]; then
  revision=${version_line#*\(}
  revision=${revision%\)}
else
  warn "the binary does not record the commit it was built from"
fi
# Unlike a development version such as 2.1.0.dev0, the git describe output
# names a commit. As in the archives GitHub makes of tags, its "v" is dropped.
package_version=${revision#v}
package_version=${package_version:-$version}
if [[ ! $package_version =~ ^[0-9A-Za-z][0-9A-Za-z.+_-]*$ ]]; then
  die "cannot name the package after: $version_line"
fi
if [[ $revision == *-dirty ]]; then
  warn "the binary was built from a checkout with uncommitted changes"
  revision=${revision%-dirty}
fi
if [[ $revision == *-g* ]]; then
  revision=${revision##*-g}
fi

src_dir=$(cache_path pono_SOURCE_DIR)
commit=$revision
if [[ -n $revision ]]; then
  commit=$(git -C "$src_dir" rev-parse -q --verify "$revision^{commit}" \
    2>/dev/null) || commit=$revision
fi
pono_url=https://github.com/stanford-centaur/pono
remote=$(git -C "$src_dir" remote get-url origin 2>/dev/null) || remote=
if [[ $remote =~ github\.com[:/](.+)$ ]]; then
  pono_url=https://github.com/${BASH_REMATCH[1]%.git}
fi

name=pono-$package_version
if [[ -e $out_dir/$name || -e $out_dir/$name.zip ]]; then
  die "$out_dir/$name or $out_dir/$name.zip already exists"
fi
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
pkg=$tmp/$name
licenses=$pkg/LICENSES
mkdir -p "$licenses"

nm -C "$binary" >"$tmp/symbols" 2>/dev/null || true
if [[ ! -s $tmp/symbols ]]; then
  die "$binary has no symbol table, package it before stripping"
fi
# A symbol that each component defines, to tell whether it is linked in.
declare -A probes=(
  [mathsat]=msat_create_env
  [cvc5]='cvc5::Solver::Solver('
  [bitwuzla]='bitwuzla::Bitwuzla::Bitwuzla('
  [btor2tools]=btor2parser_new
  [libpoly]=lp_polynomial_new
  [gmp]=__gmpz_init
  [mpfr]=mpfr_init2
)
# Components that other configurations link, whose licenses this script
# does not cover.
declare -A unsupported=(
  [Boolector]=boolector_new
  [Yices2]=yices_init
  [Z3]=Z3_mk_context
  [ic3ia]='ic3ia::'
  [CoreIR]='CoreIR::'
  [gperftools]=ProfilerStart
  [CLN]='cln::'
  [GLPK]=glp_create_prob
  [CoCoALib]='CoCoA::'
  [Normaliz]='libnormaliz::'
  [Kissat]=kissat_init
  [CryptoMiniSat]='CMSat::'
)
declare -A linked=()
for component in "${!probes[@]}"; do
  if grep -q -F -e "${probes[$component]}" "$tmp/symbols"; then
    linked[$component]=1
  fi
done
for component in "${!unsupported[@]}"; do
  if grep -q -F -e "${unsupported[$component]}" "$tmp/symbols"; then
    die "the binary links $component, whose license this script misses"
  fi
done

# Rows of name, version, license, files, and source URL, separated by the
# ASCII unit separator, which unlike a tab never merges adjacent fields.
sep=$'\x1f'
components=()

# Adds the component with fields $1 to $5 to the rows.
add_component() {
  local IFS=$sep
  components+=("$*")
}

# Records in array $1 that solver $3 builds that dependency at revision $2.
add_user() {
  local -n users=$1
  users[$2]="${users[$2]:+${users[$2]} and }$3"
}

cp "$src_dir/LICENSE" "$pkg/LICENSE"
add_component Pono "$version_line" BSD-3-Clause LICENSE \
  "$pono_url${commit:+/tree/$commit}"

if [[ -v linked[mathsat] ]]; then
  msat_root=$(cache_path MathSAT_ROOT)
  mkdir "$licenses/mathsat"
  cp "$msat_root/LICENSE.txt" "$msat_root/CREDITS.txt" "$licenses/mathsat/"
  # The library records no version; the executable beside it reports it.
  msat_version=
  if [[ -x $msat_root/bin/mathsat ]]; then
    msat_version=$("$msat_root/bin/mathsat" -version |
      sed -n '1s/^MathSAT5 version \([^,]*\),.*/\1)/p')
  fi
  add_component "MathSAT 5" "${msat_version:-unknown}" \
    "LicenseRef-MathSAT5 (research and evaluation only)" \
    "LICENSES/mathsat/LICENSE.txt, LICENSES/mathsat/CREDITS.txt" \
    https://mathsat.fbk.eu
fi

# contrib/setup-smt-switch.sh installs smt-switch into its source tree.
smt_switch_src=$(dirname "$(cache_path SmtSwitch_ROOT)")
[[ -f $smt_switch_src/LICENSE ]] || die "$smt_switch_src/LICENSE not found"
cp "$smt_switch_src/LICENSE" "$licenses/smt-switch.txt"
smt_switch_version=$(git -C "$smt_switch_src" describe --tags --always \
  --dirty 2>/dev/null) || smt_switch_version=unknown
if [[ $smt_switch_version == *-dirty ]]; then
  warn "smt-switch in $smt_switch_src has uncommitted changes"
fi
smt_switch_url=https://github.com/stanford-centaur/smt-switch
smt_switch_commit=$(git -C "$smt_switch_src" rev-parse HEAD 2>/dev/null) ||
  smt_switch_commit=
add_component smt-switch "$smt_switch_version" BSD-3-Clause \
  LICENSES/smt-switch.txt \
  "$smt_switch_url${smt_switch_commit:+/tree/$smt_switch_commit}"

# Revisions of CaDiCaL and SymFPU, each mapped to the solvers that build it.
declare -A cadical=() symfpu=()

if [[ -v linked[cvc5] ]]; then
  cvc5_dir=$(cache_path cvc5_DIR)
  cvc5_version=$(sed -n 's/^set(PACKAGE_VERSION "\(.*\)")$/\1/p' \
    "$cvc5_dir/cvc5ConfigVersion.cmake")
  [[ -n $cvc5_version ]] || die "cannot read the cvc5 version in $cvc5_dir"
  cvc5_tag=cvc5-$cvc5_version
  mkdir "$licenses/cvc5"
  github_file cvc5/cvc5 "$cvc5_tag" COPYING "$licenses/cvc5/COPYING"
  github_file cvc5/cvc5 "$cvc5_tag" licenses/minisat-LICENSE \
    "$licenses/cvc5/minisat-LICENSE"
  add_component cvc5 "$cvc5_version" "BSD-3-Clause AND MIT" \
    "LICENSES/cvc5/COPYING, LICENSES/cvc5/minisat-LICENSE" \
    "https://github.com/cvc5/cvc5/tree/$cvc5_tag"
  # cvc5 downloads the dependencies it builds at the versions that its CMake
  # find modules pin.
  for module in CaDiCaL SymFPU Poly; do
    github_file cvc5/cvc5 "$cvc5_tag" "cmake/Find$module.cmake" \
      "$tmp/Find$module.cmake"
  done
  tag=$(cmake_set CaDiCaL_VERSION "$tmp/FindCaDiCaL.cmake")
  add_user cadical "$tag" cvc5
  rev=$(cmake_set SymFPU_COMMIT "$tmp/FindSymFPU.cmake")
  add_user symfpu "$rev" cvc5
  if [[ -v linked[libpoly] ]]; then
    poly_version=$(cmake_set Poly_VERSION "$tmp/FindPoly.cmake")
    add_component LibPoly "$poly_version" LGPL-3.0-or-later \
      "LICENSES/LGPL-3.0.txt, LICENSES/GPL-3.0.txt" \
      "https://github.com/SRI-CSL/libpoly/tree/v$poly_version"
  fi
elif [[ -v linked[libpoly] ]]; then
  die "the binary links LibPoly without cvc5, which this script misses"
fi

if [[ -v linked[bitwuzla] ]]; then
  bzla_root=$(cache_path Bitwuzla_ROOT)
  # contrib/setup-smt-switch.sh builds Bitwuzla at the revision smt-switch
  # pins.
  bzla_pin=$smt_switch_src/cmake/provision/bitwuzla/pin.cmake
  [[ -f $bzla_pin ]] || die "$bzla_pin not found"
  bzla_rev=$(sed -n 's/^ *\(COMMIT\|TAG\) \(.*\)$/\2/p' "$bzla_pin")
  [[ -n $bzla_rev ]] || die "$bzla_pin pins no revision"
  bzla_version=${bzla_rev:0:12}
  bzla_pc=$(find "$bzla_root" -name bitwuzla.pc -print -quit)
  if [[ -n $bzla_pc ]]; then
    bzla_version="$(sed -n 's/^Version: //p' "$bzla_pc") ($bzla_version)"
  fi
  github_file bitwuzla/bitwuzla "$bzla_rev" COPYING "$licenses/bitwuzla.txt"
  add_component Bitwuzla "$bzla_version" MIT LICENSES/bitwuzla.txt \
    "https://github.com/bitwuzla/bitwuzla/tree/$bzla_rev"
  # Bitwuzla builds the CaDiCaL and SymFPU that its Meson wraps pin.
  for wrap in cadical symfpu; do
    github_file bitwuzla/bitwuzla "$bzla_rev" "subprojects/$wrap.wrap" \
      "$tmp/$wrap.wrap"
  done
  tag=$(sed -n 's|^source_url = .*/archive/\(.*\)\.tar\.gz$|\1|p' \
    "$tmp/cadical.wrap")
  rev=$(sed -n 's/^revision = //p' "$tmp/symfpu.wrap")
  [[ -n $tag && -n $rev ]] || die "cannot read the wraps of Bitwuzla"
  add_user cadical "$tag" Bitwuzla
  add_user symfpu "$rev" Bitwuzla
fi

for tag in "${!cadical[@]}"; do
  github_file arminbiere/cadical "$tag" LICENSE "$licenses/cadical-$tag.txt"
  add_component "CaDiCaL (in ${cadical[$tag]})" "$tag" MIT \
    "LICENSES/cadical-$tag.txt" \
    "https://github.com/arminbiere/cadical/tree/$tag"
done

for rev in "${!symfpu[@]}"; do
  dir=symfpu-${rev:0:12}
  mkdir "$licenses/$dir"
  for file in LICENSE LICENSE-BSD; do
    github_file martin-cs/symfpu "$rev" "$file" "$licenses/$dir/$file"
  done
  add_component "SymFPU (in ${symfpu[$rev]})" "${rev:0:12}" \
    "BSD-3-Clause, chosen from GPL-3.0 OR BSD-3-Clause" \
    "LICENSES/$dir/LICENSE, LICENSES/$dir/LICENSE-BSD" \
    "https://github.com/martin-cs/symfpu/tree/$rev"
done

if [[ -v linked[btor2tools] ]]; then
  # contrib/setup-btor2tools.sh clones btor2tools beside its install prefix.
  btor2_src=$(dirname "$(cache_path Btor2Tools_ROOT)")/btor2tools
  if [[ -e $btor2_src/.git ]]; then
    btor2_rev=$(git -C "$btor2_src" rev-parse HEAD)
  else
    btor2_rev=$(sed -n 's/^BTOR2TOOLS_VERSION=//p' \
      "$src_dir/contrib/setup-btor2tools.sh")
  fi
  [[ -n $btor2_rev ]] || die "cannot tell which btor2tools was built"
  github_file hwmcc/btor2tools "$btor2_rev" LICENSE.txt \
    "$licenses/btor2tools.txt"
  add_component btor2tools "${btor2_rev:0:12}" MIT LICENSES/btor2tools.txt \
    "https://github.com/hwmcc/btor2tools/tree/$btor2_rev"
fi

# Pono's CMake fetches fmt into the build directory.
fmt_src=$(cache_path FMT_SOURCE_DIR)
cp "$fmt_src/LICENSE" "$licenses/fmt.txt"
fmt_number=$(sed -n '/^#define FMT_VERSION /{s/[^0-9]*//p;q;}' \
  "$fmt_src"/include/fmt/*.h)
[[ -n $fmt_number ]] || die "cannot read the fmt version in $fmt_src"
fmt_version=$((fmt_number / 10000)).$((fmt_number / 100 % 100))
fmt_version+=.$((fmt_number % 100))
add_component fmt "$fmt_version" MIT LICENSES/fmt.txt \
  "https://github.com/fmtlib/fmt/tree/$fmt_version"

# The license of the option parser is the comment its header opens with.
optionparser=("$src_dir"/contrib/optionparser-*/src/optionparser.h)
[[ -f ${optionparser[0]} ]] || die "optionparser.h not found in $src_dir"
sed -n '1,/\*\//p' "${optionparser[0]}" >"$licenses/optionparser.txt"
optionparser_version=${optionparser[0]#"$src_dir"/contrib/optionparser-}
add_component "The Lean Mean C++ Option Parser" \
  "${optionparser_version%%/*}" MIT LICENSES/optionparser.txt \
  http://optionparser.sourceforge.net

if [[ -v linked[gmp] ]]; then
  gmp_h=$(cache_path GMP_INCLUDE_DIR)/gmp.h
  gmp_version=$(sed -n \
    's/^#define __GNU_MP_VERSION\(_MINOR\|_PATCHLEVEL\)\? *\([0-9]*\)$/\2/p' \
    "$gmp_h" | paste -s -d .)
  add_component GMP "$gmp_version" \
    "LGPL-3.0-or-later, chosen from LGPL-3.0-or-later OR GPL-2.0-or-later" \
    "LICENSES/LGPL-3.0.txt, LICENSES/GPL-3.0.txt" https://gmplib.org
fi

if [[ -v linked[mpfr] ]]; then
  mpfr_h=$(cache_path MPFR_INCLUDE_DIR)/mpfr.h
  mpfr_version=$(sed -n 's/^#define MPFR_VERSION_STRING "\(.*\)"$/\1/p' \
    "$mpfr_h")
  add_component MPFR "$mpfr_version" LGPL-3.0-or-later \
    "LICENSES/LGPL-3.0.txt, LICENSES/GPL-3.0.txt" https://www.mpfr.org
fi

# A static glibc records no version in the binary, so this assumes the
# binary was linked on this machine.
glibc_version=$(getconf GNU_LIBC_VERSION)
add_component "GNU C Library" "${glibc_version##* }" LGPL-2.1-or-later \
  LICENSES/LGPL-2.1.txt https://www.gnu.org/software/libc/

comment=$(readelf -p .comment "$binary")
gcc_version=$(sed -n 's/^.*GCC: (.*) \([0-9.]*\)$/\1/p' <<<"$comment" |
  sort -V | tail -n 1)
if [[ -n $gcc_version ]]; then
  github_file gcc-mirror/gcc "refs/tags/releases/gcc-$gcc_version" \
    COPYING.RUNTIME "$licenses/GCC-exception-3.1.txt"
  add_component "GCC runtime libraries (libstdc++, libgcc)" "$gcc_version" \
    "GPL-3.0-or-later WITH GCC-exception-3.1" \
    "LICENSES/GPL-3.0.txt, LICENSES/GCC-exception-3.1.txt" \
    https://gcc.gnu.org
fi

# Copies GNU license $1 from the system, or downloads it from path $2 on
# gnu.org, to LICENSES/$3.
gnu_license() {
  if [[ -f /usr/share/common-licenses/$1 ]]; then
    cp "/usr/share/common-licenses/$1" "$licenses/$3"
  else
    download "https://www.gnu.org/licenses/$2" "$licenses/$3"
  fi
}
gnu_license GPL-3 gpl-3.0.txt GPL-3.0.txt
gnu_license LGPL-3 lgpl-3.0.txt LGPL-3.0.txt
gnu_license LGPL-2.1 old-licenses/lgpl-2.1.txt LGPL-2.1.txt

rule=$(printf '%79s' '' | tr ' ' -)
{
  printf '\n\n%s\n Software Linked into This Binary\n%s\n\n' "$rule" "$rule"
  cat <<'EOF'
This binary statically links the software below, each under its own
license. Their license texts are in the LICENSES directory.
EOF
  if [[ -v linked[mathsat] ]]; then
    cat <<'EOF'

MathSAT 5 is available for research and evaluation purposes only. It can
not be used in a commercial environment, particularly as part of a
commercial product, without written permission from Fondazione Bruno
Kessler. This restriction therefore applies to the binary as a whole.
EOF
  fi
  # The first row is pono itself, which the text above already covers.
  for row in "${components[@]:1}"; do
    IFS=$sep read -r c_name c_version c_license c_files c_url <<<"$row"
    printf '\n%s %s\n' "$c_name" "$c_version"
    printf '  License: %s\n  Files:   %s\n  Source:  %s\n' \
      "$c_license" "$c_files" "$c_url"
  done
} >>"$pkg/LICENSE"

configure_flags=--static
if [[ $build_type == Debug ]]; then
  configure_flags+=" --debug"
fi
if [[ -v linked[mathsat] ]]; then
  configure_flags+=" --with-msat"
fi
built_from="<$pono_url>"
if [[ -n $commit ]]; then
  built_from="commit \`$commit\` of $built_from"
fi
component_table=$(
  for row in "${components[@]}"; do
    IFS=$sep read -r c_name c_version c_license c_files c_url <<<"$row"
    printf '| %s | %s | %s | <%s> |\n' \
      "$c_name" "$c_version" "$c_license" "$c_url"
  done
)
# Fill in README.md.in: keep the lines between @IF_MATHSAT@ and @END_MATHSAT@
# only when MathSAT is linked, then replace each @NAME@ placeholder. Each
# value is quoted so that bash does not read & in it as the matched text.
readme=$(<"$assets/README.md.in")
if [[ -v linked[mathsat] ]]; then
  readme=$(sed '/^@\(IF\|END\)_MATHSAT@$/d' <<<"$readme")
else
  readme=$(sed '/^@IF_MATHSAT@$/,/^@END_MATHSAT@$/d' <<<"$readme")
fi
readme=${readme//@PACKAGE_VERSION@/"$package_version"}
readme=${readme//@VERSION_LINE@/"$version_line"}
readme=${readme//@COMPONENTS@/"$component_table"}
readme=${readme//@SOURCE@/"$built_from"}
readme=${readme//@CONFIGURE_FLAGS@/"$configure_flags"}
if [[ $readme =~ @[A-Z_]+@ ]]; then
  die "README.md.in has an unknown placeholder: ${BASH_REMATCH[0]}"
fi
printf '%s\n' "$readme" >"$pkg/README.md"

if [[ $strip_binary == yes ]]; then
  strip --strip-all -o "$pkg/pono" "$binary"
else
  cp "$binary" "$pkg/pono"
fi
cp "$assets/smoketest.sh" "$pkg/"
cp -r "$assets/examples" "$pkg/"
chmod -R u=rwX,go=rX "$pkg"
# SV-COMP runs the smoke test in the package directory, but the tool has to
# work from any directory, so run it from another one.
if ! (cd "$tmp" && "$pkg/smoketest.sh") >"$tmp/smoketest.log" 2>&1; then
  cat "$tmp/smoketest.log" >&2
  die "the smoke test of the package fails"
fi

mkdir -p "$out_dir"
mv "$pkg" "$out_dir/"
(cd "$out_dir" && zip -q -r -X "$name.zip" "$name")
echo "Wrote $out_dir/$name and $out_dir/$name.zip"
if [[ -n $commit ]]; then
  echo "Check that $pono_url has commit $commit, which README.md cites."
fi
