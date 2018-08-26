.PHONY: all clean

all: sample sample-weighted sample-non_weighted sample-fast

STD_CXXFLAGS = -std=c++11
#STD_CXXFLAGS = -std=c++14

CXXFLAGS += $(STD_CXXFLAGS)

newton_method/libnewton.a:
	$(MAKE) -C newton_method

%: %.cc

%: %.o newton_method/libnewton.a
	$(LINK.cc) $^ $(LOADLIBES) $(LDLIBS) -o $@

clean:
	$(MAKE) -C newton_method clean
	-rm -f *~ *.o
