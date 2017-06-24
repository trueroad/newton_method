.PHONY: all clean

all: sample sample-weighted sample-non_weighted

STD_CXXFLAGS = -std=c++11
#STD_CXXFLAGS = -std=c++14

EIGEN_CXXFLAGS = $(shell pkg-config eigen3 --cflags)
#EIGEN_CXXFLAGS = -I./eigen33

#DEBUG_CXXFLAGS =
DEBUG_CXXFLAGS = -DDEBUG_NEWTON_METHOD

CXXFLAGS = $(STD_CXXFLAGS) $(EIGEN_CXXFLAGS) $(DEBUG_CXXFLAGS)

sample: sample.o newton.o
	$(LINK.cc) $^ $(LOADLIBES) $(LDLIBS) -o $@

sample-weighted: sample-weighted.o newton.o
	$(LINK.cc) $^ $(LOADLIBES) $(LDLIBS) -o $@

sample-non_weighted: sample-non_weighted.o newton.o
	$(LINK.cc) $^ $(LOADLIBES) $(LDLIBS) -o $@

clean:
	-rm -f *~ *.o
