# set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 600, 400 
# set output 'surface1.9.png'
set samples 21, 21
set isosamples 11, 11
set view 130,20,1
set style data lines
set title "3D gnuplot demo" 
set xlabel "X axis" 
set xlabel  offset character -3, -2, 0 font "" textcolor lt -1 norotate
set xrange [ 0.0000 : 1.5700 ] noreverse nowriteback
set ylabel "Y axis" 
set ylabel  offset character 3, -2, 0 font "" textcolor lt -1 rotate
set yrange [ 0.0000 : 1.04700 ] noreverse nowriteback
set zlabel "Z axis" 
set zlabel  offset character -5, 0, 0 font "" textcolor lt -1 norotate
set zrange [ 0.5 : 1.70 ] noreverse nowriteback
DEBUG_TERM_HTIC = 119
DEBUG_TERM_VTIC = 119
## Last datafile plotted: "$grid"
splot "wave.dat" u 2:3:4 with lines,  \
#  splot "wave.dat" u 2:3:4 with points pt 1 ps -1,  \
#"wave.dat" u 2:3:5 with points,  \

   
