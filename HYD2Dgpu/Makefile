##########################################
# Makefile for Turbulence Studies
##########################################

##########################################
# Programs
##########################################

exe=Simulation.x
ana=Analysis.x

######################	
# complile options
######################
fc= nvfortran
foptopenacc = -Minfo=accel -acc
fopt = -g -traceback -O2
#fopt = -g -traceback -check all -fpe0

######################
# directory
######################

dird := bindata
dira := output

all: ${exe} ${ana}

.PHONY: all clean allclean

#################
# Fortran
#################


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
	rm -fr ${dird} ${dira}
