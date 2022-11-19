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
 
set style line 91 lt 1 lw 6 lc rgb "black" # 
set style line 92 lt 2 lw 6 lc rgb "black" #

# input file
ifnum=600
input= sprintf("output/spc%05d.dat",ifnum)

##########################################
# Evolusiton of shock and Gain Radius
##########################################

if(pngflag==1)set output "k-E_k.png"
set log 

set xlabel "Wave number"
set xrange [1:100]

set ylabel "E_k"
set yrange [*:*]

set key right top

plot NaN notitle \
, input  u 1:2  notitle w l ls 1  \
, 6*x**(-5.0/3.0) title "-5/3" w l ls 2


if(pngflag==1)set output "k-V_k.png"
set log 

set ylabel "Enstrophy"

plot NaN notitle \
, input  u 1:3  notitle w l ls 1  \

reset
set term pop
