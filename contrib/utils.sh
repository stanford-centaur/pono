# shellcheck shell=bash
# Helpers shared by the dependency setup scripts. Source it, do not run it.

# Aborts unless file $1 has SHA256 digest $2, removing path $3 when given.
verify_sha256() {
  local actual
  if command -v sha256sum >/dev/null 2>&1; then
    actual="$(sha256sum "$1" | cut -d ' ' -f 1)"
  else
    # macOS has no sha256sum, but ships shasum with its Perl.
    actual="$(shasum -a 256 "$1" | cut -d ' ' -f 1)"
  fi
  if [[ $actual != "$2" ]]; then
    echo "Checksum mismatch for $1" >&2
    echo "  expected SHA256: $2" >&2
    echo "  actual SHA256:   $actual" >&2
    if [[ -n ${3-} ]]; then
      rm -rf "$3"
    fi
    exit 1
  fi
}
