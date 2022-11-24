##########################################
# Makefile for Turbulence Studies
##########################################

##########################################
# Programs
##########################################

exe=Simulation.x
ana=Analysis.x
vis=Visualization.x

######################	
# complile options
######################
fc=ifort -extend-source
fopt=-g -traceback -O2
#fopt=-g -traceback -check all -fpe0

##########################################
# Timeseqencial data and figures
##########################################

######################
# directory
######################

dird := bindata
dira := output
dirf := figures
dirm := movies

countfile=control.dat

######################
# files
######################

filename = ${dird}/unf00001.dat
fileexists = $(shell ls | grep ${filename})
ifeq (${fileexists}, ${filename})
    BIN= $(shell ls ${dird}/unf*.dat)
else
    BIN= ${filename}
endif

VOR   = $(patsubst ${dird}/unf%.dat,${dira}/vor%.dat,$(BIN))
VORPNG= $(patsubst ${dira}/vor%.dat,${dirf}/vor%.png,$(VOR))
KINPNG= $(patsubst ${dira}/vor%.dat,${dirf}/kin%.png,$(VOR))

SPC   = $(patsubst ${dird}/unf%.dat,${dira}/spc%.dat,$(BIN))
KSPPNG= $(patsubst ${dira}/spc%.dat,${dirf}/ksp%.png,$(SPC))
HKSPNG= $(patsubst ${dira}/spc%.dat,${dirf}/hks%.png,$(SPC))

##########################################
# Movies
##########################################

kinmovie=${dirm}/anikin.mp4
vormovie=${dirm}/anivor.mp4
kspmovie=${dirm}/aniksp.mp4
hksmovie=${dirm}/anihks.mp4

moviefile=${vormovie} ${kinmovie} ${kspmovie} ${hksmovie}

all: ${exe} ${ana} ${vis} ${moviefile}

.PHONY: all clean allclean

#################
# ffmpeg
#################

movie: ${moviefile}

${vormovie}: MakeMovie.sh ${VORPNG}
	./MakeMovie.sh vor

${kinmovie}: MakeMovie.sh ${KINPNG}
	./MakeMovie.sh kin

${kspmovie}: MakeMovie.sh ${KSPPNG}
	./MakeMovie.sh ksp

${hksmovie}: MakeMovie.sh ${HKSPNG}


#################
# gnuplot
#################

spectrum: ${specfile}

${specfile}: spectrum.plt
	gnuplot spectrum.plt

2Dsnaps: ${VORPNG}

${VORPNG}: vor.plt MakeSnap.sh ${VOR}
	./MakeSnap.sh vor.plt

#################
# visualization
#################

HDF5 = /home/takiwkkz/hdf5
FLAGHDF =  -I${HDF5}/include -L${HDF5}/lib -lhdf5 -lhdf5_fortran 

${vis}: Visualization.f90
	${fc} ${fopt} ${FLAGHDF} $< -o ${vis}

${VOR}: ${vis} ${BIN} ${countfile}
	./${vis}

#################
# analysis
#################

analysis-code: ${ana} 

count-number: ${countfile}

run-analysis: ${ana} ${count}
	./${ana}

${countfile}: Init.sh
	./Init.sh

${SPC}: ${ana} ${BIN} ${countfile}
	./${ana}


#################
# simulation
#################

run: ${BIN}
	./${exe}

${BIN}: ${exe}
	./${exe}


simulation-code: ${exe} 

${ana}: Analysis.f90
	${fc} ${fopt} $< -o ${ana}

${exe}: Simulation.f90
	${fc} ${fopt} $< -o ${exe}

#################
# clean up
#################
clean:
	rm -f ${exe} ${ana} *.o *.mod *~

allclean:
	rm -fr ${dirm} ${dirf} ${dird} ${dira} ${count}


