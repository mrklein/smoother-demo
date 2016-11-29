SHELL = /bin/bash
SOURCES = $(wildcard *.cc)
OBJECTS = $(SOURCES:.cc=.o)

EXECUTABLE = ../bin/smoother

FOAM_HOME = $(HOME)/OpenFOAM/OpenFOAM-4.x

# Local definitions
DEFS = darwin64 \
       WM_ARCH_OPTION=64 \
       WM_DP WM_LABEL_SIZE=32 \
       NoRepository
$(foreach def,$(DEFS),$(eval CPPFLAGS += -D$(def)))

# These are always in CPPFLAGS
CPPFLAGS += -I$(FOAM_HOME)/src/OpenFOAM/lnInclude \
	    -I$(FOAM_HOME)/src/OSspecific/POSIX/lnInclude
# Libraries location
LDFLAGS += -L$(FOAM_HOME)/platforms/darwin64ClangDPInt32Opt/lib

MODULES = finiteVolume meshTools
$(foreach m,$(MODULES),$(eval CPPFLAGS += -I$(FOAM_HOME)/src/$(m)/lnInclude))
$(foreach m,$(MODULES),$(eval LDFLAGS += -l$(m)))

CXXFLAGS += -std=c++14

LDFLAGS += -Wl,-execute,-undefined,dynamic_lookup \
	   -lOpenFOAM

.SUFFIXES: .cc .o
.PHONY: clean
.DEFAULT: all

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $^

clean:
	$(RM) $(OBJECTS) $(EXECUTABLE)

