dird := bindata
dira := output
dirf := figures
dirm := movies

MOVIES := ${dirm}/anivor.mp4

BIN= ${dird}/unf00001.dat $(shell ls ${dird}/unf*.dat)
VOR= $(patsubst ${dird}/unf%.dat,${dira}/vor%.dat,$(BIN))
VORPNG=$(patsubst ${dira}/vor%.dat,${dirf}/vor%.png,$(VOR))

exe=Simulation.x
ana=Analysis.x

all: ${exe} ${MOVIES}

.PHONY: all clean allclean

#################
# ffmpeg
#################
${dirm}/anivor.mp4: MakeMovie.sh ${VORPNG}
	./MakeMovie.sh vor

#################
# gnuplot
#################

${VORPNG}: vor.plt MakeSnap.sh ${VOR}
	rm -fr ${dirf}/vor*.png
	./MakeSnap.sh vor.plt

#################
# analysis
#################

${VOR}: ${ana} ${BIN} control.dat
	rm -fr ${dira}/vor*.dat
	./${ana}

control.dat: Init.sh ${BIN}
	./Init.sh

#################
# simulation
#################

${BIN}: ${exe}
	rm -fr ${dird}/unf*.dat
	rm -fr ${dird}/bin*.dat
	./${exe}


fc=ifort -extend-source
fopt=-g -traceback -O2
#fopt=-g -traceback -check all -fpe0

${ana}: dataanalysis.f90
	${fc} ${fopt} $< -o ${ana}

${exe}: main.f90
	${fc} ${fopt} $< -o ${exe}

#################
# clean up
#################
clean:
	rm -f ${exe} ${ana} *.o *.mod *~

allclean:
	rm -fr ${dirm} ${dirf} ${dird} ${dirf}
