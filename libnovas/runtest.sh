#!/bin/sh

OUTFILE="testnovas.out"
export LD_LIBRARY_PATH=`pwd`
./testnovas > $OUTFILE

diff -u Cdist/checkout-stars-usno.txt $OUTFILE
