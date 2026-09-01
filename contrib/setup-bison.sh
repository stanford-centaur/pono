#!/usr/bin/env bash
set -euo pipefail

version=3.8.2
url=https://mirrors.ocf.berkeley.edu/gnu/bison
dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd)/deps/bison
if [[ -d $dir ]]; then
  echo "$dir already exists; remove it if you want to redownload bison" && exit
fi
archive=bison-$version.tar.gz
curl -LsS -o $archive $url/$archive
mkdir -p "$dir"
tar -xf $archive -C "$dir" --strip-components 1
rm -f $archive
cd "$dir" && ./configure --prefix "$dir-install" && make && make install
