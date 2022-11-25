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

# 2D cross section
VOR   = $(patsubst ${dird}/unf%.dat,${dira}/vor%.dat,$(BIN))
VORPNG= $(patsubst ${dira}/vor%.dat,${dirf}/vor%.png,$(VOR))
KINPNG= $(patsubst ${dira}/vor%.dat,${dirf}/kin%.png,$(VOR))

# spectrum
SPC   = $(patsubst ${dird}/unf%.dat,${dira}/spc%.dat,$(BIN))
KSPPNG= $(patsubst ${dira}/spc%.dat,${dirf}/ksp%.png,$(SPC))
VSPPNG= $(patsubst ${dira}/spc%.dat,${dirf}/vsp%.png,$(SPC))

# intgrated values 
TOT   = $(patsubst ${dird}/unf%.dat,${dira}/tot%.dat,$(BIN))

##########################################
# Movies
##########################################

kinmovie=${dirm}/anikin.mp4
vormovie=${dirm}/anivor.mp4

2Dmovies=${vormovie} ${kinmovie}

kspmovie=${dirm}/aniksp.mp4
vspmovie=${dirm}/anivsp.mp4

specmovies= ${kspmovie} ${vspmovie}


moviefile=${2Dmovies} ${specmovies}
##########################################
# Final spectrum
##########################################

specfile=k-E_k.png

##########################################
# Time evolution
##########################################

timefile=t-E.png

##########################################
# Proceadures
##########################################

all: ${exe} ${moviefile} ${specfile} ${timefile}

.PHONY: all clean allclean

##########################################
# Time evolution
##########################################
${timefile}: TimeEvolution.plt t-prof.dat
	gnuplot TimeEvolution.plt

t-prof.dat: MakeTimeseq.sh ${TOT}
	./MakeTimeseq.sh

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


########################
# 2D cross section
########################

2Dsnaps: ${VORPNG} ${KINPNG} 

${VORPNG}: vor.plt MakeSnap.sh ${VOR}
	./MakeSnap.sh vor.plt

${KINPNG}: vor.plt MakeSnap.sh ${VOR}
	./MakeSnap.sh vor.plt

########################
# spectrum
########################

spectrum: ${specfile}

${specfile}: spectrum-fin.plt
	gnuplot spectrum-fin.plt

specs: ${KSPPNG} ${VSPPNG} 

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
