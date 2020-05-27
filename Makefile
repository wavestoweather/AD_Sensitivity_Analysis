GCC=g++
PCC=scorep-g++
GCCFLAGS= --std=c++14 -pthread -lgsl -lgslcblas -lm -DCODI_UseForcedInlines -fargument-noalias-global -ftree-loop-vectorize -lnetcdf_c++4 -lnetcdf #-DLIKWID_PERFMON #  -ffast-math
GCCINCLUDES=-I.
TIMESTEPPER=-DRK4ICE
ATMOFLAGS=-DCONSTANT_DROP=FALSE
SEASON=-DSPRING
FLUX=-DFLUX
SOURCE=-DWCB2 -DSAT_CALC

BUILD=build
OBJ_DIR=$(BUILD)/objects
APP_DIR=$(BUILD)/apps

SRC=\
    $(wildcard src/microphysics/*.cpp) \

SRC_SCRATCH=src/scratch/load_test.cpp

SRC_SCAN=src/scratch/scan.cpp

OBJECTS=$(SRC:%.cpp=$(OBJ_DIR)/%.o)
TARGETS=$(SRC:%.cpp=$(APP_DIR)/%)

OBJECTS_SCRATCH=$(SRC_SCRATCH:%.cpp=$(OBJ_DIR)/%.o)
TARGETS_SCRATCH=$(SRC_SCRATCH:%.cpp=$(APP_DIR)/%)

OBJECTS_SCAN=$(SRC_SCAN:%.cpp=$(OBJ_DIR)/%.o)
TARGETS_SCAN=$(SRC_SCAN:%.cpp=$(APP_DIR)/%)

all: build $(TARGETS)

scratch: build $(TARGETS_SCRATCH)

scan: build $(TARGETS_SCAN)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(GCC) $(GCCINCLUDES) $(GCCFLAGS) $(TIMESTEPPER) $(SEASON) $(FLUX) $(SOURCE) -o $@ -c $<

$(TARGETS): $(OBJECTS)
	@mkdir -p $(@D)
	$(GCC) $(GCCINCLUDES) $(GCCFLAGS) -o $@ $<

$(TARGETS_SCRATCH): $(OBJECTS_SCRATCH)
	@mkdir -p $(@D)
	$(GCC) $(GCCINCLUDES) $(GCCFLAGS) -o $@ $<

$(TARGETS_SCAN): $(OBJECTS_SCAN)
	@mkdir -p $(@D)
	$(GCC) $(GCCINCLUDES) $(GCCFLAGS) -o $@ $<


.PHONY: all build clean debug release

build:
	@mkdir -p $(BUILD)
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(APP_DIR)

debug: GCCFLAGS += -g -Og
debug: all

debug_scratch: GCCFLAGS += -g -Og
debug_scratch: scratch

release: GCCLFAGS += -O2 -march=native
release: all

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/src/microphysics/*

