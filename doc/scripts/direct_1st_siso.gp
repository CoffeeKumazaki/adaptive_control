set term pngcairo dashed font "Courier,24" size 1280, 960 ``

path = system("dirname ".ARG0)."/"
set output path."../img/direct_1st_siso.png"
load path."linetype.gp"

set xlabel "t"
set ylabel "y"
p path."../../build/direct_1st_siso.txt" u 1:2 w l ti "y_{plant}", \
  "" u 1:3 dt 2 w l ti "y_{model}", \
  "" u 1:($3-$2) w l ti "error" lt 2 lw 1 lc "black"

set output

set term pngcairo dashed font "Courier,24" size 1280, 960 ``

set output path."../img/direct_1st_siso_prm.png"

set xlabel "t"
set ylabel "parameters"
p path."../../build/direct_1st_siso.txt" u 1:5 w l ti "k_y", \
  "" u 1:6 dt 1 w l ti "k_r", \
  -1.5 lt 1 dt 2 ti "", \
  1 lt 2 dt 2 ti ""

set output