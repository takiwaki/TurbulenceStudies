reset

pngflag=1
if(pngflag==1)set terminal push
if(pngflag==1)set terminal pngcairo dashed enhanced font "Helvetica,20"
if(pngflag==1)set encoding utf8

set style line 1 lt 1 lw 6 lc rgb "#ff2800" # universal design red 
set style line 2 lt 1 lw 6 lc rgb "#0041ff" # universal design blue
set style line 3 lt 1 lw 6 lc rgb "#35a16B" # universal design green
set style line 4 lt 1 lw 6 lc rgb "#faf500" # universal design yellow
set style line 5 lt 1 lw 6 lc rgb "#66ccff" # universal design sky-blue,azure
set style line 6 lt 1 lw 6 lc rgb "#ff99a0" # universal design pink
set style line 7 lt 1 lw 6 lc rgb "#ff9900" # universal design orange
set style line 8 lt 1 lw 6 lc rgb "#9a0079" # universal design purple
set style line 9 lt 1 lw 6 lc rgb "#663300" # universal design brown
 
set style line 91 lt 1 lw 2 lc rgb "black" # 
set style line 92 lt 2 lw 2 lc rgb "black" #

# input file
ifnum=580
inputHYD2D= sprintf("../HYD2Dgpu/output/spc%05d.dat",ifnum)
inputMHD2D= sprintf("../MHD2Dgpu/output/spc%05d.dat",ifnum)
inputHYD2Dlow= sprintf("../HYD2D/output/spc%05d.dat",ifnum)
inputMHD2Dlow= sprintf("../MHD2D/output/spc%05d.dat",ifnum)

ifnum=581
inputHYD3D= sprintf("../HYD3Dgpu/output/spc%05d.dat",ifnum)
inputMHD3D= sprintf("../MHD3Dgpu/output/spc%05d.dat",ifnum)
inputHYD3Dlow= sprintf("../HYD2D/output/spc%05d.dat",ifnum)
inputMHD3Dlow= sprintf("../MHD2D/output/spc%05d.dat",ifnum)


if(pngflag==1)set output "k-EK_k2D.png"
set log 

set xlabel "Wave number"
set xrange [1:100]

set ylabel "Spectrum [a.u.]"
set format y "10^{%L}"
set yrange [1.0e-3:1.0e1]

set key right top

plot NaN notitle \
, inputHYD2D  u 1:2      title "Kinetic energy (2D  HYD)" w l ls 1  \
, inputMHD2D  u 1:($2/5) title "Kinetic energy (2D  MHD)" w l ls 7  \
, inputMHD2D  u 1:($6/10)   title "Magnetic energy (2D  MHD)" w l ls 2  \
, 0.1*(x/20.0)**(-5.0/3.0)   title "-5/3"   w l ls 91 \
#, 0.1*(x/10.0)**(-1.0)   title "-1"   w l ls 92


if(pngflag==1)set output "k-EK_k3D.png"
set log 

set xlabel "Wave number"
set xrange [1:100]

set yrange [*:*]

set key right top

plot NaN notitle \
, inputHYD3D  u 1:2    title "Kinetic energy (3D  HYD)" w l ls 1  \
, inputMHD3D  u 1:($2) title "Kinetic energy (3D  MHD)" w l ls 7  \
, inputMHD3D  u 1:($6*10) title "Magnetic energy (3D  MHD)" w l ls 2  \
, 4.0*(x/10.0)**(-5.0/3.0)   title "-5/3"   w l ls 91 \
#, 0.1*(x/10.0)**(-1.0)   title "-1"   w l ls 92





reset
set term pop


