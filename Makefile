GCC=g++
PCC=scorep-g++
GCCFLAGS= --std=c++14 -pthread -lgsl -lgslcblas -lm -DCODI_UseForcedInlines -fargument-noalias-global -ftree-loop-vectorize -lnetcdf_c++4 -lnetcdf #-DLIKWID_PERFMON #  -ffast-math
GCCINCLUDES=-L/mnt/localscratch/lib/
TIMESTEPPER=-DRK4ICE
ATMOFLAGS=-DCONSTANT_DROP=FALSE
SEASON=-DSPRING
FLUX=-DFLUX
SOURCE=-DWCB

all: trajectories score debug test trajectories_sb trajectories_sb_noice

score: trajectories.cpp
	$(PCC) $(GCCINCLUDES) trajectories.cpp -o trajectories_score.out  $(GCCFLAGS) $(TIMESTEPPER) $(SEASON) $(FLUX) $(SOURCE)

trajectories: trajectories.cpp
	$(GCC) $(GCCINCLUDES) trajectories.cpp -o trajectories.out  $(GCCFLAGS) -O2 -march=native -DRK4 $(SEASON) $(FLUX) $(SOURCE)

trajectories_sb: trajectories.cpp
	$(GCC) $(GCCINCLUDES) trajectories.cpp -o trajectories_sb.out  $(GCCFLAGS) -O2 -march=native -DRK4ICE $(SEASON) $(FLUX) $(SOURCE)

trajectories_sb_noice: trajectories.cpp
	$(GCC) $(GCCINCLUDES) trajectories.cpp -o trajectories_sb_noice.out  $(GCCFLAGS) -O2 -march=native -DRK4NOICE $(SEASON) $(FLUX) $(SOURCE)

debug: trajectories.cpp
	$(GCC) $(GCCINCLUDES) trajectories.cpp -o trajectories_debug.out  $(GCCFLAGS) -g -Og $(TIMESTEPPER) $(SEASON) $(FLUX) $(SOURCE)

test: load_test.cpp
	$(GCC) $(GCCINCLUDES) load_test.cpp -o load_test.out  $(GCCFLAGS) -O2 -march=native $(TIMESTEPPER) $(SEASON) $(FLUX) $(SOURCE)

clean:
	rm -f *.out