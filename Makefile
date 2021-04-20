GCC=c++
PCC=scorep-g++
#
GCCFLAGS=--std=c++14 -lgsl -lgslcblas -lm -DCODI_UseForcedInlines -fargument-noalias-global -ftree-loop-vectorize -lnetcdf_c++4 -lnetcdf -Wall #-DLIKWID_PERFMON #  -ffast-math
GCCINCLUDES=-I.
TIMESTEPPER=-DRK4ICE
ATMOFLAGS=-DCONSTANT_DROP=FALSE
SEASON=-DAUTUMN
FLUX=-DFLUX
SOURCE=-DMET3D -DSB_CONV -DSB_SHAPE -DNPROCS=4

BUILD=build
OBJ_DIR=$(BUILD)/objects
APP_DIR=$(BUILD)/apps

SRC=\
    $(wildcard src/microphysics/*.cpp) \

SRC_SCRATCH=src/scratch/load_test.cpp

SRC_SCAN=src/scratch/scan.cpp

SRC_NETCDF=src/scratch/netcdf_test.cpp

OBJECTS=$(SRC:%.cpp=$(OBJ_DIR)/%.o)
TARGETS=$(SRC:%.cpp=$(APP_DIR)/%)

OBJECTS_SCRATCH=$(SRC_SCRATCH:%.cpp=$(OBJ_DIR)/%.o)
TARGETS_SCRATCH=$(SRC_SCRATCH:%.cpp=$(APP_DIR)/%)

OBJECTS_SCAN=$(SRC_SCAN:%.cpp=$(OBJ_DIR)/%.o)
TARGETS_SCAN=$(SRC_SCAN:%.cpp=$(APP_DIR)/%)

OBJECTS_NETCDF=$(SRC_NETCDF:%.cpp=$(OBJ_DIR)/%.o)
TARGETS_NETCDF=$(SRC_NETCDF:%.cpp=$(APP_DIR)/%)

all: build $(TARGETS)

scratch: build $(TARGETS_SCRATCH)

scan: build $(TARGETS_SCAN)

netcdf: build $(TARGETS_NETCDF)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(GCC) $(GCCINCLUDES) $(SOURCE) $(GCCFLAGS) $(TIMESTEPPER) $(SEASON) $(FLUX) -o $@ -c $<

$(TARGETS): $(OBJECTS)
	@mkdir -p $(@D)
	$(GCC) $(GCCINCLUDES) $(GCCFLAGS) -o $@ $<

$(TARGETS_SCRATCH): $(OBJECTS_SCRATCH)
	@mkdir -p $(@D)
	$(GCC) $(GCCINCLUDES) $(GCCFLAGS) -o $@ $<

$(TARGETS_SCAN): $(OBJECTS_SCAN)
	@mkdir -p $(@D)
	$(GCC) $(GCCINCLUDES) $(GCCFLAGS) -o $@ $<

$(TARGETS_NETCDF): $(OBJECTS_NETCDF)
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

debug_netcdf: GCCFLAGS += -g -Og
debug_netcdf: netcdf

clean:
	-@trash-put $(OBJ_DIR)/*
	-@trash-put $(APP_DIR)/src/microphysics/*

