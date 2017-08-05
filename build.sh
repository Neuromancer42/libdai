#!/usr/bin/env bash

if [[ $HOSTNAME =~ "fir" ]]; then
   make -j 8 CCINC="-I $HOME/Error-Ranking/chord-fork/libs/boost_1_61_0 -I $HOME/Error-Ranking/chord-fork/libs/gmp-6.1.0 -I include " CCLIB="-L $HOME/Error-Ranking/chord-fork/libs/boost_1_61_0/installdir/lib -L $HOME/Error-Ranking/chord-fork/libs/gmp-6.1.0/installdir/lib -L lib " HOSTNAME=$HOSTNAME
   g++ -std=c++11 -O2 -march=native -Wall -Wextra -Werror -I ./include -I $HOME/Error-Ranking/chord-fork/libs/boost_1_61_0 -I $HOME/Error-Ranking/chord-fork/libs/gmp-6.1.0 -L $HOME/Error-Ranking/chord-fork/libs/gmp-6.1.0/installdir/lib -lgmp -lgmpxx wrapper.cpp ./lib/libdai.a -o wrapper
else
   make -j 8 && \
   g++ -std=c++11 -O2 -march=native -Wall -Wextra -Werror -I ./include -lgmp -lgmpxx wrapper.cpp ./lib/libdai.a -o wrapper
fi

exit $?
