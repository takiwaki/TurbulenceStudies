dird := bindata
dira := output
dirf := figures
dirm := movies

MOVIES := ${dirm}/anidnx.mp4

BIN= ${dird}/unf00001.dat $(shell ls ${dird}/unf*.dat)
VOR= $(patsubst ${dird}/unf%.dat,${dira}/vor%.dat,$(BIN))

exe=Simulation.x
ana=Analysis.x

all: ${exe} ${MOVIES}

.PHONY: all clean allclean

#################
# ffmpeg
#################
${dirm}/anidnx.mp4: MakeMovie.sh ${DNPNG} Init.sh
	./Init.sh
	./MakeMovie.sh dnx

#################
# gnuplot
#################
${VORPNG}: vor.plt MakeSnap.sh ${VOR}
	rm -fr ${dirf}/vor*.png
	./MakeSnap.sh vor.plt

${VOR}: ${ana} ${BIN} control.dat
	rm -fr ${dira}/vor*.dat
	./${ana}


#################
# simulation
#################

control.dat: Init.sh ${BIN}
	./Init.sh

${BIN}: ${exe}
	rm -fr ${dird}/unf*.dat
	rm -fr ${dird}/bin*.dat
	./${exe}


#################
# Fortran
#################

fc=ifort -extend-source
fopt=-g -traceback -O2
#fopt=-g -traceback -check all -fpe0

${exe}: main.F
	${fc} ${fopt} $< -o ${exe}

${ana}: dataanalysis.f90
	${fc} ${fopt} $< -o ${ana}

#################
# clean up
#################
clean:
	rm -f ${exe} *.o *.mod *~

allclean:
	rm -fr ${dirm} ${dirf} ${dird} ${dirf}
