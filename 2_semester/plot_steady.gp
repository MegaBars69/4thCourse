# Плотность
set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'density_steady.png'
set view map
set xlabel 'x'
set ylabel 'y'
set title 'Плотность газа {/Symbol r}'
set size ratio -1
set xrange [0:3]
set yrange [0:2]
set palette defined (0 "blue", 1 "red")
splot 'result_steady.dat' using 1:2:3 with pm3d notitle

# Поле скорости
set output 'velocity_steady.png'
set title 'Поле скорости'
unset view
plot 'result_steady.dat' using 1:2:4:5 with vectors head filled lc rgb 'black' notitle
