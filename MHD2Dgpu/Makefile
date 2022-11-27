dird := bindata

exe=Simulation.x
ana=Analysis.x

all: ${exe} ${ana}

.PHONY: all clean allclean

#################
# Fortran
#################

fc= nvfortran
foptopenacc = -Minfo=accel -acc
fopt = -g -traceback -O2
#fopt = -g -traceback -check all -fpe0

${exe}: Simulation.f90
	${fc} ${fopt} ${foptopenacc} $< -o ${exe}

${ana}: Analysis.f90
	${fc} ${fopt} ${foptopenacc} $< -o ${ana}

#################
# clean up
#################
clean:
	rm -f ${exe} *.o *.mod *~

allclean:
	rm -fr ${dird}
