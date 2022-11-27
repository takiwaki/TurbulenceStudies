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
#fopt=-g -traceback -check all -fpe0

HDF5 = /home/takiwkkz/hdf5
FLAGHDF =  -I${HDF5}/include -L${HDF5}/lib -lhdf5 -lhdf5_fortran 

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
JCDPNG= $(patsubst ${dira}/vor%.dat,${dirf}/jcd%.png,$(VOR))
MAGPNG= $(patsubst ${dira}/vor%.dat,${dirf}/mag%.png,$(VOR))

SPC   = $(patsubst ${dird}/unf%.dat,${dira}/spc%.dat,$(BIN))
HKSPNG= $(patsubst ${dira}/spc%.dat,${dirf}/hks%.png,$(SPC))
KSPPNG= $(patsubst ${dira}/spc%.dat,${dirf}/ksp%.png,$(SPC))
HCRPNG= $(patsubst ${dira}/spc%.dat,${dirf}/hcr%.png,$(SPC))
HMMPNG= $(patsubst ${dira}/spc%.dat,${dirf}/hmm%.png,$(SPC))
MSPPNG= $(patsubst ${dira}/spc%.dat,${dirf}/msp%.png,$(SPC))

##########################################
# Movies
##########################################

# 2D snaps evolution
vormovie=${dirm}/anivor.mp4
kinmovie=${dirm}/anikin.mp4
jcdmovie=${dirm}/anijcd.mp4
magmovie=${dirm}/animag.mp4

# spectrum evolution
kspmovie=${dirm}/aniksp.mp4
hksmovie=${dirm}/anihks.mp4
hcrmovie=${dirm}/anihcr.mp4
hmmmovie=${dirm}/anihmm.mp4
mspmovie=${dirm}/animsp.mp4

moviefile=${vormovie} ${kinmovie} ${jcdmovie} ${magmovie} ${kspmovie} ${hksmovie} ${hcrmovie} ${hmmmovie} ${mspmovie}

all: ${exe} ${moviefile}

.PHONY: all clean allclean

#################
# ffmpeg
#################

movie: ${moviefile}

${vormovie}: MakeMovie.sh ${VORPNG}
	./MakeMovie.sh vor

${kinmovie}: MakeMovie.sh ${KINPNG}
	./MakeMovie.sh kin

${jcdmovie}: MakeMovie.sh ${JCDPNG}
	./MakeMovie.sh jcd

${magmovie}: MakeMovie.sh ${MAGPNG}
	./MakeMovie.sh mag

# 5 kinds of spectrum

${kspmovie}: MakeMovie.sh ${KSPPNG}
	./MakeMovie.sh ksp

${hksmovie}: MakeMovie.sh ${HKSPNG}
	./MakeMovie.sh hks

${hcrmovie}: MakeMovie.sh ${HCRPNG}
	./MakeMovie.sh hcr

${hmmmovie}: MakeMovie.sh ${MPTPNG}
	./MakeMovie.sh hmm

${mspmovie}: MakeMovie.sh ${MSPPNG}
	./MakeMovie.sh msp

#################
# gnuplot
#################

spectrum: ${specfile}

${specfile}: spectrum.plt
	gnuplot spectrum.plt

2Dsnaps: ${VORPNG}

${VORPNG}: vor.plt MakeSnap.sh ${VOR}
	./MakeSnap.sh vor.plt

${JCDPNG}: vor.plt MakeSnap.sh ${VOR}
	./MakeSnap.sh vor.plt

${KINPNG}: vor.plt MakeSnap.sh ${VOR}
	./MakeSnap.sh vor.plt

${MAGPNG}: vor.plt MakeSnap.sh ${VOR}
	./MakeSnap.sh vor.plt


${KSPPNG}: vor.plt MakeSnap.sh ${SPC}
	./MakeSnap.sh spectrum.plt

${HKSPNG}: spectrum.plt MakeSnap.sh ${SPC}
	./MakeSnap.sh spectrum.plt

${HCRPNG}: spectrum.plt MakeSnap.sh ${SPC}
	./MakeSnap.sh spectrum.plt

${MSPPNG}: spectrum.plt MakeSnap.sh ${SPC}
	./MakeSnap.sh spectrum.plt

${HMMPNG}: spectrum.plt MakeSnap.sh ${SPC}
	./MakeSnap.sh spectrum.plt


#################
# visualization
#################


${vis}: Visualization.f90
	${fc} ${fopt} ${FLAGHDF} $< -o ${vis}

${VOR}: ${vis} ${BIN} ${countfile}
	./${vis}

#################
# analysis
#################

count-number: ${countfile}

${countfile}: CountBindata.sh
	./CountBindata.sh

${SPC}: ${ana} ${BIN} ${countfile}
	./${ana}

#################
# simulation
#################

${BIN}: ${exe}
	./${exe}

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


