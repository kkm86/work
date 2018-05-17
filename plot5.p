# Gnuplot script file for plotting data in file "res.dat"
      # This file is called   plot5.p
 set termopt enhanced    # turn on enhanced text mode
   
      
   
      set   autoscale                        # scale axes automatically
      unset log                              # remove any log-scaling
      unset label                            # remove any previous labels
      unset arrow
      unset xtic
      set title font "Times,20"
      set xlabel font "Times,15"
      set ylabel font "Times,15"	 
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically
      set key font ",12"
      set title "Radial electron probability density in Potassium 3d/4s"
      set xlabel "Bohr radius"
      set ylabel "4{/Symbol p}r^2{/Symbol r}(r)"
	 set xr [0.0:4.0]
	 set yr [0.0:25.1]
	 plot    "res3.dat" using 2:3 title'3d'  with lines,  \
	         "res1.dat" using 2:3 title'4s'  with lines,  \
# "res3.dat" using 2:3 title'3'  with linespoint,  \
	         
