dird := bindata
dira := output
dirf := figures
dirm := movies

BIN= ${dird}/unf00001.dat $(shell ls ${dird}/unf*.dat)
VOR= $(patsubst ${dird}/unf%.dat,${dira}/vor%.dat,$(BIN))
SPC= $(patsubst ${dird}/unf%.dat,${dira}/spc%.dat,$(BIN))
VORPNG=$(patsubst ${dira}/vor%.dat,${dirf}/vor%.png,$(VOR))

exe=Simulation.x
ana=Analysis.x
vis=Visualization.x

countfile=control.dat
moviefile=${dirm}/anivor.mp4
specfile=k-E_k.png

#################
# complile options
#################
fc=ifort -extend-source
fopt=-g -traceback -O2
#fopt=-g -traceback -check all -fpe0

all: ${exe} ${moviefile}

.PHONY: all clean allclean

#################
# ffmpeg
#################

movie: ${moviefile}

${moviefile}: MakeMovie.sh ${VORPNG}
	./MakeMovie.sh vor

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


