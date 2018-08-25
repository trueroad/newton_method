#!/bin/sh

REPLACE_SV=`echo $1 | sed -e "s/\([^\-]\+\)-\([^\-]\+\)-\([^\-]\+\)\.cc/\1/"`
REPLACE_LS=`echo $1 | sed -e "s/\([^\-]\+\)-\([^\-]\+\)-\([^\-]\+\)\.cc/\2/"`
REPLACE_AL=`echo $1 | sed -e "s/\([^\-]\+\)-\([^\-]\+\)-\([^\-]\+\)\.cc/\3/"`

cat ${REPLACE_SV}-LEAST_SQUARE-ALGORITHM.cc.in | \
    sed -e "s/LEAST_SQUARE/${REPLACE_LS}/g" | \
    sed -e "s/ALGORITHM/${REPLACE_AL}/g" > $1
