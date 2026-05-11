# Плотность
set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'density.png'
set view map
set xlabel 'x'
set ylabel 'y'
set title 'Плотность газа {/Symbol r}'
set size ratio -1
set pm3d map
splot 'result.dat' using 1:2:3 with pm3d notitle

# Поле скорости
set output 'velocity.png'
set title 'Поле скорости'
plot 'result.dat' using 1:2:4:5 with vectors head filled lc rgb 'black' notitle