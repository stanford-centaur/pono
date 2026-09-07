#!/usr/bin/env bash
set -e

dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
dir="$dir/../deps/mathsat"

usage() {
  cat <<EOF
Usage: $0 [<option> ...]

Downloads the MathSAT SMT Solver. Note this solver is under a custom (non BSD compliant) license.

-h, --help              display this message and exit
-y, --auto-yes          automatically agree to conditions (default: off)
EOF
  exit 0
}

die() {
  echo "*** $0: $*" >&2
  exit 1
}

# Aborts unless archive $1 has the expected SHA256 digest $2, discarding
# $dir so that a rerun starts over. macOS has no sha256sum, so fall back
# to the shasum that ships with its Perl.
verify_sha256() {
  local actual
  if command -v sha256sum >/dev/null 2>&1; then
    actual="$(sha256sum "$1" | cut -d ' ' -f 1)"
  else
    actual="$(shasum -a 256 "$1" | cut -d ' ' -f 1)"
  fi
  if [[ $actual != "$2" ]]; then
    echo "Checksum mismatch for $msat_archive"
    echo "  expected SHA256: $2"
    echo "  actual SHA256:   $actual"
    rm -rf "$dir"
    exit 1
  fi
}

get_msat=default
msat_version="5.6.12"
msat_release_url="https://mathsat.fbk.eu/release/mathsat"
# Digests of the $msat_version archives. Update them together with
# msat_version, otherwise the download stops verifying.
msat_sha256_linux_x86_64="\
1de984ed8500ce0895970116d57f73b381929fa6600d4e29b6beb3042f2b721b"
msat_sha256_macos="\
28fe0711bdd920217af706b07983f033470efc2ea15a4563397e1a930b75e9f5"

while [[ $# -gt 0 ]]; do
  case $1 in
    -h | --help) usage ;;
    -y | --auto-yes) get_msat=y ;;
    *) die "unexpected argument: $1" ;;
  esac
  shift
done

if [[ $get_msat == default ]]; then
  read -rp "MathSAT is distributed under a custom (non-BSD compliant) ""\
license. By continuing you acknowledge this distinction and assume ""\
responsibility for meeting the license conditions. Continue? ""\
[y]es/[n]o: " get_msat
fi

if [[ $get_msat != y ]]; then
  echo "Not downloading MathSAT"
  exit 0
fi

if [[ ! -d $dir ]]; then
  # Each platform picks its archive and that archive's digest together.
  if [[ $OSTYPE == linux* || $OSTYPE == cygwin* ]]; then
    msat_archive=$msat_release_url-$msat_version-linux-x86_64.tar.gz
    msat_sha256=$msat_sha256_linux_x86_64
  elif [[ $OSTYPE == darwin* ]]; then
    msat_archive=$msat_release_url-$msat_version-macos.tar.gz
    msat_sha256=$msat_sha256_macos
  else
    echo "Unrecognized OSTYPE=$OSTYPE"
    exit 1
  fi

  # Catches a platform case added without a digest to pin its archive.
  if [[ -z $msat_sha256 ]]; then
    echo "No SHA256 digest is pinned for $msat_archive"
    exit 1
  fi

  mkdir -p "$dir"
  cd "$dir"
  wget -O mathsat.tar.gz "$msat_archive"
  verify_sha256 mathsat.tar.gz "$msat_sha256"

  tar -xf mathsat.tar.gz --strip-components 1
  rm mathsat.tar.gz

else
  echo "$dir already exists. If you want to re-download, please remove it."
fi

if [[ -f "$dir/lib/libmathsat.a" ]]; then
  echo "It appears mathsat was setup successfully into $dir."
  echo "You may now install it with ./configure.sh --msat && cd build && make"
else
  echo "Downloading mathsat failed."
  echo "Please see their website: http://mathsat.fbk.eu/"
  exit 1
fi
