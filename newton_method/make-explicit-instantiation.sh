#!/bin/sh

REPLACE_LS=`echo $1 | sed -e "s/solve-\([^\-]\+\)-\([^\-]\+\)\.cc/\1/"`
REPLACE_AL=`echo $1 | sed -e "s/solve-\([^\-]\+\)-\([^\-]\+\)\.cc/\2/"`

cat solve-LEAST_SQUARE-ALGORITHM.cc.in | \
    sed -e "s/LEAST_SQUARE/${REPLACE_LS}/g" | \
    sed -e "s/ALGORITHM/${REPLACE_AL}/g" > $1
