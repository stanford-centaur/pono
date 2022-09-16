#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DEPS=$DIR/../deps

usage () {
cat <<EOF
Usage: $0 [<option> ...]

Configures the CMAKE build environment.

-h, --help              display this message and exit
--msat-home             mathsat location (default: deps/mathsat)
EOF
  exit 0
}

die () {
    echo "*** $0: $*" 1>&2
    exit 1
}

mkdir -p $DEPS
cd $DEPS

msat_home="$DEPS/mathsat"
ic3ia_version=ic3ia-22.07

while [ $# -gt 0 ]
do
    case $1 in
        -h|--help) usage;;
        --msat-home=*)
            msat_home=${1##*=}
            # Check if msat_home is an absolute path and if not, make it
            # absolute.
            case $msat_home in
                /*) ;;                            # absolute path
                *) msat_home=$(pwd)/$msat_home ;; # make absolute path
            esac
            ;;
        *) die "unexpected argument: $1";;
    esac
    shift
done

if [ -d "$DEPS/ic3ia" ]; then
    echo "It appears ic3ia was already downloaded"
    echo "If you'd like to rebuild, please manually delete it"
    exit 1
fi

if [ ! -d "$msat_home" ]; then
    echo "Could not locate mathsat"
    echo "Required for building ic3ia"
    echo "Looked here: $msat_home"
    exit 1
fi

curl -Lk https://es-static.fbk.eu/people/griggio/ic3ia/$ic3ia_version.tar.gz --output $DEPS/$ic3ia_version.tar.gz

if [ ! -f "$DEPS/$ic3ia_version.tar.gz" ]; then
    echo "It appears that downloading ic3ia failed."
    exit 1
fi

tar -xf $ic3ia_version.tar.gz
rm $ic3ia_version.tar.gz
mv $ic3ia_version ic3ia
cd ic3ia
mkdir build
cd build
cmake .. -DMATHSAT_DIR=$msat_home -DCMAKE_BUILD_TYPE=Release
make -j
cd $DEPS/..
