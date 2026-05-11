# Корректное совмещение плотности и векторов скорости
set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set xlabel 'x'
set ylabel 'y'
set size ratio -1
set xrange [0:3]
set yrange [0:2]
set view map
set palette defined (0 "blue", 1 "red")

set output 'frame_T1.png'
set title 'Момент времени T/4'
splot 'result_T2.50.dat' using 1:2:3 with pm3d notitle, \
      'result_T2.50.dat' using 1:2:(0.1):4:5:(0.0) with vectors head filled lc rgb 'black' notitle

set output 'frame_T2.png'
set title 'Момент времени T/2'
splot 'result_T5.00.dat' using 1:2:3 with pm3d notitle, \
      'result_T5.00.dat' using 1:2:(0.1):4:5:(0.0) with vectors head filled lc rgb 'black' notitle

set output 'frame_T3.png'
set title 'Момент времени T/3*4'
splot 'result_T7.50.dat' using 1:2:3 with pm3d notitle, \
      'result_T7.50.dat' using 1:2:(0.1):4:5:(0.0) with vectors head filled lc rgb 'black' notitle

set output 'frame_T4.png'
set title 'Момент времени T/1'
splot 'result_T10.00.dat' using 1:2:3 with pm3d notitle, \
      'result_T10.00.dat' using 1:2:(0.1):4:5:(0.0) with vectors head filled lc rgb 'black' notitle

