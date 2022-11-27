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
set style line 92 lt 2 lw 6 lc rgb "black" #

# input file
if (exist("ifnum")==0 ) ifnum=581
input= sprintf("output/spc%05d.dat",ifnum)

##########################################
# Kinetic energy
##########################################

outputfile="k-E_k.png"
#outputfile= sprintf("figures/ksp%05d.png",ifnum)
if(pngflag==1)set output outputfile

set log 
set format y "10^{%L}"

set xlabel "Wave number"
set xrange [1:100]

set ylabel "Kinetic energy"
set yrange [*:*]

set key right top

plot NaN notitle \
, input  u 1:2  notitle w l ls 1  \
, (0.2)*(x/10)**(-5.0/3.0) title "-5/3" w l ls 91

##########################################
# Kinetic helicity
##########################################

outputfile="k-HK_k.png"
#outputfile= sprintf("figures/hsk%05d.png",ifnum)
if(pngflag==1)set output outputfile
set log 
set format y "10^{%L}"

set ylabel "Kinetic helicity"

plot NaN notitle \
, input  u 1:(abs($3))  notitle w l ls 1  \

##########################################
# Both
##########################################

outputfile="k-E_kHK_k.png"
if(pngflag==1)set output outputfile
set ylabel "Spectrum [a.u.]"

set log 
set format y "10^{%L}"

plot NaN notitle \
, input  u 1:(5*$2)    title "Kinetic energy"   w l ls 1  \
, input  u 1:(abs($3)) title "Kinetic helicity" w l ls 2  \
, (20.0)*(x/10)**(-5.0/3.0) title "-5/3" w l ls 91

reset
set term pop


