all: libnewton.a

#
# Options
#
#CPPFLAGS_NEWTON_METHOD =
CPPFLAGS_NEWTON_METHOD = -DDEBUG_NEWTON_METHOD


#
# Objects and libraries
#
OBJS_COMMON = newton.o
OBJS_SOLVE = \
	solve-through_pass-ColPivHouseholderQR.o \
	solve-weighted-ColPivHouseholderQR.o \

LDLIBS =
CPPFLAGS_EIGEN = $(shell pkg-config eigen3 --cflags)
#CPPFLAGS_EIGEN = -I./eigen33


#
# Programs and other settings
#
#RANLIB = ranlib
CXXFLAGS_STD = -std=c++11
#CXXFLAGS_STD = -std=c++14


#
# Dependencies
#
DEPS = $(OBJS:.o=.d)
CPPFLAGS += -MMD -MP -MF $(@:.o=.d) -MT $@
-include $(DEPS)


#
# Flags
#
CXXFLAGS += $(CXXFLAGS_STD)
CPPFLAGS += $(CPPFLAGS_EIGEN) $(CPPFLAGS_NEWTON_METHOD)

OBJS = $(OBJS_COMMON) $(OBJS_SOLVE)
SOURCES_SOLVE = $(OBJS_SOLVE:.o=.cc)


#
# Rules
#

clean:
	$(RM) *~ $(OBJS) $(DEPS) $(SOURCES_SOLVE)

.PHONY: all clean

libnewton.a: $(OBJS)
	$(RM) $@
	$(AR) $(ARFLAGS)c $@ $^
#	$(RANLIB) $@

$(SOURCES_SOLVE): solve-LEAST_SQUARE-ALGORITHM.cc.in
	./make-explicit-instantiation.sh $@

#
# Implicit variables (from `make -p`)
#

#CXX = g++
#AR = ar
#ARFLAGS = rv
#RM = rm -f
#COMPILE.cc = $(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c
#OUTPUT_OPTION = -o $@


#
# Implicit rules (from `make -p`)
#

#%.o: %.cc
#	$(COMPILE.cc) $(OUTPUT_OPTION) $<
