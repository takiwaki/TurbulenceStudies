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
JCDPNG= $(patsubst ${dira}/vor%.dat,${dirf}/jcd%.png,$(VOR))
MAGPNG= $(patsubst ${dira}/vor%.dat,${dirf}/mag%.png,$(VOR))

SPC   = $(patsubst ${dird}/unf%.dat,${dira}/spc%.dat,$(BIN))
VSPPNG= $(patsubst ${dira}/spc%.dat,${dirf}/vsp%.png,$(SPC))
KSPPNG= $(patsubst ${dira}/spc%.dat,${dirf}/ksp%.png,$(SPC))
HCRPNG= $(patsubst ${dira}/spc%.dat,${dirf}/hcr%.png,$(SPC))
MPTPNG= $(patsubst ${dira}/spc%.dat,${dirf}/mpt%.png,$(SPC))
MSPPNG= $(patsubst ${dira}/spc%.dat,${dirf}/msp%.png,$(SPC))

##########################################
# Movies
##########################################

vormovie=${dirm}/anivor.mp4
kinmovie=${dirm}/anikin.mp4

jcdmovie=${dirm}/anijcd.mp4
magmovie=${dirm}/animag.mp4

kspmovie=${dirm}/aniksp.mp4
vspmovie=${dirm}/anivsp.mp4
hcrmovie=${dirm}/anihcr.mp4
mptmovie=${dirm}/animpt.mp4
mspmovie=${dirm}/animsp.mp4

moviefile=${vormovie} ${kinmovie} ${jcdmovie} ${magmovie} ${kspmovie} ${vspmovie} ${hcrmovie} ${mptmovie} ${mspmovie}

specfile=k-E_k.png

all: ${exe} ${moviefile} ${specfile}

.PHONY: all clean allclean

#################
# ffmpeg
#################

movie: ${jcdmovie} ${vormovie} ${kinmovie} ${magmovie} 

${vormovie}: MakeMovie.sh ${VORPNG}
	./MakeMovie.sh vor

${jcdmovie}: MakeMovie.sh ${JCDPNG}
	./MakeMovie.sh jcd

${kinmovie}: MakeMovie.sh ${VORPNG}
	./MakeMovie.sh kin

${magmovie}: MakeMovie.sh ${MAGPNG}
	./MakeMovie.sh mag


# 5 kinds of spectrum

${kspmovie}: MakeMovie.sh ${KSPPNG}
	./MakeMovie.sh ksp

${vspmovie}: MakeMovie.sh ${VSPPNG}
	./MakeMovie.sh vsp

${hcrmovie}: MakeMovie.sh ${HCRPNG}
	./MakeMovie.sh hcr

${mptmovie}: MakeMovie.sh ${MPTPNG}
	./MakeMovie.sh mpt

${mspmovie}: MakeMovie.sh ${MSPPNG}
	./MakeMovie.sh msp


#################
# gnuplot
#################

spectrum: ${specfile}

${specfile}: spectrum.plt
	gnuplot spectrum.plt

2Dsnaps: ${VORPNG} ${JCDPNG} ${KINPNG} ${MAGPNG}

${VORPNG}: vor.plt MakeSnap.sh ${VOR}
	./MakeSnap.sh vor.plt

${JCDPNG}: vor.plt MakeSnap.sh ${VOR}
	./MakeSnap.sh vor.plt

${KINPNG}: vor.plt MakeSnap.sh ${VOR}
	./MakeSnap.sh vor.plt

${MAGPNG}: vor.plt MakeSnap.sh ${VOR}
	./MakeSnap.sh vor.plt


${KSPPNG}: spectrum.plt MakeSnap.sh ${SPC}
	./MakeSnap.sh spectrum.plt

${VSPPNG}: spectrum.plt MakeSnap.sh ${SPC}
	./MakeSnap.sh spectrum.plt

${HCRPNG}: spectrum.plt MakeSnap.sh ${SPC}
	./MakeSnap.sh spectrum.plt

${MSPPNG}: spectrum.plt MakeSnap.sh ${SPC}
	./MakeSnap.sh spectrum.plt

${MPTPNG}: spectrum.plt MakeSnap.sh ${SPC}
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
