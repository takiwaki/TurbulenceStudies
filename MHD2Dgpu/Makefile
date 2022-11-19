dird := bindata
dira := output
dirf := figures
dirm := movies

MOVIES := ${dirm}/anivor.mp4 ${dirm}/anijcd.mp4

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
	./MakeMovie.sh jcd

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
	rm -fr ${dird} ${dira} ${dirf} ${dirm} control.dat
