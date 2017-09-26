#!/usr/bin/env bash
curdir=`pwd`
if [[ $HOSTNAME =~ "fir" ]]; then
   # fir machines do not have boost and libgmp installed. Therefore, we use pre-compiled libraries that we install in chord-fork/libs.
   # Also, not all the libraries of boost are available in our installation (ex: -lboost_unit_test_framework is not available). 
   # Therefore, we suppress the build of unittests, utils, etc and build only the libdai library, if the build happens on a fir machine. This check is done using HOSTNAME, in the Makefile 
   make -j 8 CCINC="-I ../libs/boost_1_61_0 -I ../libs/gmp-6.1.0 -I include " CCLIB="-L ../libs/boost_1_61_0/installdir/lib -L ../libs/gmp-6.1.0/installdir/lib -L lib " HOSTNAME=$HOSTNAME
   g++ -std=c++11 -O2 -march=native -Wall -Wextra -Werror -I ./include -I ../libs/boost_1_61_0 -I ../libs/gmp-6.1.0 -Wl,-rpath -Wl,$curdir/../libs/gmp-6.1.0/installdir/lib -L ../libs/gmp-6.1.0/installdir/lib -lgmp -lgmpxx -fopenmp wrapper.cpp ./lib/libdai.a -o wrapper
else
   make -j 8 && \
   g++ -std=c++11 -O2 -march=native -Wall -Wextra -Werror -I ./include -lgmp -lgmpxx -fopenmp wrapper.cpp ./lib/libdai.a -o wrapper
fi

exit $?
