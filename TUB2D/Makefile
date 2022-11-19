dird := bindata
dira := output
dirf := figures
dirm := movies

MOVIES := ${dirm}/anivor.mp4

BIN= ${dird}/unf00001.dat $(shell ls ${dird}/unf*.dat)
VOR= $(patsubst ${dird}/unf%.dat,${dira}/vor%.dat,$(BIN))
VORPNG=$(patsubst ${dira}/vor%.dat,${dirf}/vor%.png,$(VOR))

exe=simplehyd.x
ana=analysis.x

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

${VOR}: ${ana} ${BIN} control.dat
	rm -fr ${dira}/vor*.dat
	./${ana}

control.dat: Init.sh ${BIN}
	./Init.sh

${BIN}: ${exe}
	rm -fr ${dird}/unf*.dat
	rm -fr ${dird}/bin*.dat
	./${exe}


fc=ifort -extend-source
fopt=-g -traceback -O2
#fopt=-g -traceback -check all -fpe0

${ana}: dataanalysis.f90
	${fc} ${fopt} $< -o ${ana}

${exe}: main.F
	${fc} ${fopt} $< -o ${exe}

#################
# clean up
#################
clean:
	rm -f ${exe} *.o *.mod *~

allclean:
	rm -f ${dirm}/*.wmv ${dirf}/*.png ${dird}/*.xss
