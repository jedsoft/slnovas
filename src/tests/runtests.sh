#!/bin/sh

export LD_LIBRARY_PATH=`pwd`/../../libnovas
OUTFILE=test.out
GOODFILE=test.ok

# valgrind --tool=memcheck --leak-check=yes --error-limit=no \
slsh ./test.sl > $OUTFILE || exit 1

diff $GOODFILE $OUTFILE
# small differences are expected.
