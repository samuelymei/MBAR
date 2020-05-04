.SUFFIXES: 
.SUFFIXES: .f90 .o

FC = ifort
FFLAGS = -CB
INCLUDE = 
LIBS = 

EXE = MBAR.x

MODULES = constant_m.mod precision_m.mod io_m.mod snapshot_m.mod simulation_m.mod reducedhamiltonian_m.mod mbar_m.mod

OBJS = precision_m.o constant.o lib.o io.o snapshot.o simulation.o reducedHamiltonian.o MBAR.o MBAR_caller.o

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

