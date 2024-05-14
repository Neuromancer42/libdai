#!/usr/bin/env bash

if [ -z "${LIBDAI}" ]; then
  echo "run `export LIBDAI=<path-to-LIBDAI-root>` before"
  exit 1
fi

if [ -z "${BINGO_RESULTS}" ]; then
  export BINGO_RESULTS="/tmp/results/"
fi
echo "set result path to ${BINGO_RESULTS}"

if [ -z $1 ]; then
  echo "usage: init.sh <real-path-to-PJBENCH>"
  exit 1
fi

if [ -L $LIBDAI/BINGO/pjbench ]; then
  echo "replace old bench link"
  rm $LIBDAI/BINGO/pjbench
fi
ln -s $1 $LIBDAI/BINGO/pjbench

# { prune-cons
pushd $LIBDAI/BINGO/scripts/bnet/prune-cons
./build.sh
if [ $? -ne 0 ]; then
  echo "Build failed: prune-cons"
  exit 1
fi
popd
# }

# { LibDAI
pushd $LIBDAI
./build.sh
if [ $? -ne 0 ]; then
  echo "Build failed: LibDAI"
  exit 1
fi
popd
# }

export PATH=$LIBDAI/BINGO/scripts/ae-wrappers/:$PATH