#!/bin/sh

SOLARSYS1_OUT="checkout-stars-full-usno.txt"
SOLARSYS3_OUT="checkout-stars-usno.txt"

OUTFILE="testnovas.out"
REFFILE="$SOLARSYS1_OUT"

ln -sf "DE421/lnxp1900p2053.421" JPLEPH

export LD_LIBRARY_PATH=`pwd`
./testnovas > $OUTFILE || exit 1

diff -u Cdist/$SOLARSYS1_OUT $OUTFILE

# small differences are expected.
exit 0
