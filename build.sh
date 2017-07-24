#!/usr/bin/env bash

make -j 8 && \
g++ -std=c++11 -O2 -march=native -Wall -Wextra -Werror -I ./include -lgmp -lgmpxx wrapper.cpp ./lib/libdai.a -o wrapper
# g++ -std=c++11 -g -Wall -Wextra -Werror -I ./include -lgmp -lgmpxx wrapper.cpp ./lib/libdai.a -o wrapper

exit $?
