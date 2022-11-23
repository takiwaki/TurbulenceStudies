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
VSPPNG= $(patsubst ${dira}/spc%.dat,${dirf}/vsp%.png,$(SPC))


##########################################
# Movies
##########################################

kinmovie=${dirm}/anikin.mp4
vormovie=${dirm}/anivor.mp4
kspmovie=${dirm}/aniksp.mp4
vspmovie=${dirm}/anivsp.mp4

moviefile=${vormovie} ${kinmovie} ${kspmovie} ${vspmovie}

specfile=k-E_k.png

all: ${exe} ${moviefile} ${specfile}

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

${vspmovie}: MakeMovie.sh ${VSPPNG}
	./MakeMovie.sh vsp

#################
# gnuplot
#################

spectrum: ${specfile}

${specfile}: spectrum.plt
	gnuplot -e ifnum=580 spectrum.plt

2Dsnaps: ${VORPNG} ${KINPNG} 

${VORPNG}: vor.plt MakeSnap.sh ${VOR}
	./MakeSnap.sh vor.plt

${KINPNG}: vor.plt MakeSnap.sh ${VOR}
	./MakeSnap.sh vor.plt

${KSPPNG}: spectrum.plt MakeSnap.sh ${SPC}
	./MakeSnap.sh spectrum.plt

${VSPPNG}: spectrum.plt MakeSnap.sh ${SPC}
	./MakeSnap.sh spectrum.plt

#################
# analysis
#################

count-number: ${countfile}

${countfile}: CountBindata.sh
	./CountBindata.sh

${VOR}: ${ana} ${BIN} ${countfile}
	./${ana}

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
	rm -fr ${dirm} ${dirf} ${dird} ${dira} ${countfile}
