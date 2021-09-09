set term pngcairo dashed font "Courier,24" size 1280, 960 ``

path = system("dirname ".ARG0)."/"
set output path."../img/direct_2nd_siso.png"
load path."linetype.gp"

set xlabel "t"
set ylabel "y"
p path."../../build/direct_2nd_siso.txt" u 1:2 w l ti "{plant}", \
  "" u 1:4 lt 2 w l ti "{model}", \
  "" u 1:($4-$2) w l ti "error" lt 2 lw 1 lc "black" 

set output

set output path."../img/direct_2nd_siso_dot.png"
load path."linetype.gp"

set xlabel "t"
set ylabel "dy/dt"
set yrange [-1:1]
p path."../../build/direct_2nd_siso.txt" u 1:3 lt 1 w l ti "{plant}", \
  "" u 1:5 lt 2 w l ti "{model}", \
  "" u 1:($5-$3) w l ti "error" lt 2 lw 1 lc "black"

set output

set term pngcairo dashed font "Courier,24" size 1280, 960 ``

set output path."../img/direct_2nd_siso_prm.png"

set xlabel "t"
set ylabel "parameters"
set yrange [*:*]
p path."../../build/direct_2nd_siso.txt" u 1:7 w l ti "k_y", \
  "" u 1:8 dt 1 w l ti "k_{dy}", \
  "" u 1:9 dt 1 w l ti "k_r", \
  0.75 lt 1 dt 2 ti "", \
  0.25 lt 2 dt 2 ti "", \
  0.25 lt 3 dt 2 ti ""

set output