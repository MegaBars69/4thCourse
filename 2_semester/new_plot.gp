# Для плотности:
set terminal pngcairo size 800,600
set output 'density_steady.png'
set view map
set size ratio -1
set xlabel 'x'
set ylabel 'y'
set title 'Плотность газа {/Symbol r}'
splot 'result_steady.dat' using 1:2:3 with pm3d notitle

# Для поля скорости:
set output 'velocity_steady.png'
set title 'Поле скорости'
plot 'result_steady.dat' using 1:2:4:5 with vectors head filled lc rgb 'black' notitle
