#####################
# initialize
######################

set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0

set size ratio -1
set view map
unset key

##########################################
# parameters
##########################################

srange=0.5

##########################################
# mainloop
##########################################

# PNG
if (exist("ifnum")==0 ) set term push
set term pngcairo enhanced font "Helvetica, 12" size 600,600
# crop 

if (exist("ifnum")==0 ) ifnum=100

ifnames = sprintf("output/vor%05d.dat",ifnum)

command = sprintf("ls %s 1> /dev/null 2> /dev/null ; echo $? ",ifnames)
flag=0
flag=system(command)
print  ifnames." found"

command = sprintf(" head -n 1 %s | sed 's/#//' ",ifnames)
time   = system(command)
print "time= ".time

# Position of color bar
set colorbox horizontal user origin 0.175, 0.92 size 0.71, 0.04
#set cbtics 1.0
#set cbtics offset 0,3.2

set size 1.0

set origin 0.05,0.0
set xlabel "X" offset 0,0
set xtics 0.25
set ylabel "Y" offset 0,0
set ytics 0.25

set xrange [-0.5:0.5]
set yrange [-0.5:0.5]


##########################
# Vorticity
##########################

set title "Vorticity"

#maxcb=20.0
#set cbrange [-maxcb:maxcb]
set cbrange [*:*]

ofname = sprintf("figures/vor%05d.png",ifnum)
set output ofname


set label 1 time at screen 0.65, screen 0.845

set palette define (-1.0 "blue", 0.0 "white", 1.0 "red")

splot  \
  ifnames u ($1):($2):3 w pm3d  \

##########################
# Current density
##########################

set title "Current density"

ofname = sprintf("figures/jcd%05d.png",ifnum)
set output ofname

#maxcb=20.0#
#set cbrange [-maxcb:maxcb]

set palette define (-1.0 "purple", 0.0 "white", 1.0 "orange")

splot  \
  ifnames u ($1):($2):4 w pm3d  \


##########################
# Kinetic energy
##########################

set title "Kinetic energy"

ofname = sprintf("figures/kin%05d.png",ifnum)
set output ofname

maxcb=2.0
#set cbrange [0:maxcb]

set palette define (0.0 "white", 1.0 "black")

splot  \
  ifnames u ($1):($2):5 w pm3d  \

##########################
# Magnetic energy
##########################

set title "Magnetic energy"

ofname = sprintf("figures/mag%05d.png",ifnum)
set output ofname

set palette define (0.0 "white", 1.0 "black")

maxcb=2.0
#set cbrange [0:maxcb]
splot  \
  ifnames u ($1):($2):6 w pm3d  \

unset label 1

reset
set term pop
