#!/usr/bin/env bash

source $(dirname "$0")/init.sh

#################################################################################
# Main Program

if [ "$#" -lt 1 ] ; then
  echo "Usage: $0 <space-separated list of benchmarks>"  >&2
  exit 1
fi

benchmarks=(${@:1})

run_benchmark() {
   bmk=$1
   if [ ! -n "${bpath[$bmk]+isset}" ]; then
      echo "Unknown benchmark: $bmk"
      exit 1
   fi
   anal=${analysis[$bmk]}
   if [ ! -d $HOME/Projects/bingo-artifact/results/$anal/$bmk/output ]; then
      mkdir -p $HOME/Projects/bingo-artifact/results/$anal/$bmk/output
   fi

   pushd $CFORK
   source ./scripts/setpaths.sh

   build-bnet.sh opt $bmk
   exact_out=$HOME/Projects/bingo-artifact/results/$anal/$bmk/output/exact
   if [ ! -d $exact_out ]; then
      mkdir -p $exact_out
      echo "Starting Bingo (opt) run..."
      ./scripts/bnet/accmd ${bnet_in[$bmk]} noaugment ${base_q[$bmk]} ${oracle_q[$bmk]} ${numiter[$bmk]} $exact_out/ &
   else echo "Bingo (opt) data is present. Not re-running. Please delete $exact_out to re-run."
   fi
   wait
   echo "$bmk done."
   popd
}

for bmk in "${benchmarks[@]}"
do
   run_benchmark $bmk &
done

wait
echo "All done. Quit."