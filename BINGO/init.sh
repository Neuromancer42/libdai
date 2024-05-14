#!/usr/bin/env bash

if [[ -z "${LIBDAI}" ]]; then
    echo "run `export LIBDAI=<path-to-LIBDAI-root>` before"
    exit 1
fi
export CFORK=$LIBDAI/BINGO

if [ -z $1 ]; then
  echo "usage: init.sh <rea-path-to-PJBENCH>"
  exit 1
fi
if [ -L $CFORK/pjbench ]; then
  echo "replace old bench link"
  rm $CFORK/pjbench
fi
ln -s $1 $CFORK/pjbench

# { prune-cons
pushd $CFORK/scripts/bnet/prune-cons
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

export PATH=$CFORK/scripts/ae-wrappers/:$PATH