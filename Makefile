.SUFFIXES: 
.SUFFIXES: .f90 .o

FC = ifort
FFLAGS = -O2 -CB
INCLUDE = 
MKL_HOME=/home/ymei/softwares/intel/mkl
LIBS = -L$(MKL_HOME)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm

EXE = MBAR.x

MODULES = constant_m.mod precision_m.mod random_m.mod io_m.mod bin_m.mod snapshot_m.mod simulation_m.mod reducedhamiltonian_m.mod mbar_m.mod

OBJS = precision_m.o constant.o lib.o random.o io.o bin.o snapshot.o simulation.o MBAR.o reducedHamiltonian.o MBAR_runner.o 

all:	${EXE}


$(EXE):$(OBJS) ${MODULES}
	$(FC) -o $@ $(FFLAGS) $(OBJS) $(LIBS)

%.o %.mod:%.f90
	$(FC) -c $(FFLAGS) $(INCLUDE) $<

include .depend

depend .depend:
	makedepf90 *.f90 > .depend

clean:
	/bin/rm -f $(EXE) $(OBJS) ${MODULES}

