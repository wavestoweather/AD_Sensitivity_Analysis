GCC=g++
PCC=scorep-g++
GCCFLAGS= --std=c++14 -pthread -lgsl -lgslcblas -lm -DCODI_UseForcedInlines -fargument-noalias-global -ftree-loop-vectorize -lnetcdf_c++4 -lnetcdf #-DLIKWID_PERFMON #  -ffast-math
GCCINCLUDES=-L/mnt/localscratch/lib/
TIMESTEPPER=-DRK4ICE
ATMOFLAGS=-DCONSTANT_DROP=FALSE
SEASON=-DSPRING
FLUX=-DFLUX
SOURCE=-DWCB

BUILD=build
OBJ_DIR=$(BUILD)/objects
APP_DIR=$(BUILD)/apps

SRC=\
    $(wildcard src/microphysics/*.cpp) \

OBJECTS=$(SRC:%.cpp=$(OBJ_DIR)/%.o)
TARGETS=$(SRC:%.cpp=$(APP_DIR)/%)

all: build $(TARGETS)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(GCC) $(GCCINCLUDES) $(GCCFLAGS) $(TIMESTEPPER) $(SEASON) $(FLUX) $(SOURCE) -o $@ -c $<

$(TARGETS): $(OBJECTS)
	@mkdir -p $(@D)
	$(GCC) $(GCCINCLUDES) $(GCCFLAGS) -o $@ $<

.PHONY: all build clean debug release

build:
	@mkdir -p $(BUILD)
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(APP_DIR)

debug: GCCFLAGS += -g -0g
debug: all

release: GCCLFAGS += -O2 -march=native
release: all

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/src/microphysics/*

