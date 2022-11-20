dird := bindata
dira := output
dirf := figures
dirm := movies

exe=Simulation.x
ana=Analysis.x
count=control.dat

#################
# complile options
#################
fc=ifort -extend-source
fopt=-g -traceback -O2
#fopt=-g -traceback -check all -fpe0

all: movie spectrum

.PHONY: all clean allclean

#################
# ffmpeg
#################
movie: MakeMovie.sh 2Dsnaps
	./MakeMovie.sh vor

#################
# gnuplot
#################

spectrum: spectrum.plt
	gnuplot spectrum.plt

2Dsnaps: vor.plt MakeSnap.sh run-analysis
	./MakeSnap.sh vor.plt

#################
# analysis
#################

analysis-code: ${ana} 

count-number: ${count}

run-analysis: ${ana} ${count} run
	./${ana}

${count}: Init.sh
	./Init.sh

#################
# simulation
#################

run: ${exe}
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
