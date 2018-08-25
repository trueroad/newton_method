# This `solve.mk` builds only the algorithms needed in samples.
# If all algorithms are required, remove this `solve.mk` then invoke `make`.
# `make-solve-ml.pl` will regenerate `solve.mk` for all algorithms.

OBJS_SOLVE = \
	solve-through_pass-ColPivHouseholderQR.o \
	solve-weighted-ColPivHouseholderQR.o \
	solve_fast-weighted-ColPivHouseholderQR.o \

