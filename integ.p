# Gnuplot script file for plotting data in file "wave.dat"
      # This file is called   integ.p
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
      set title "P_{12} and P_{21}"
	 set xlabel "log(R/a_{sc}(a.u.))"
	 set ylabel ""
	 set log x
	 set xr [0.0043:43.8]
	 set yr [-0.03:0.03]
	 plot "wave.dat" using 2:3 title 'P_{12}' with lines,  \
	    "wave.dat" using 2:4 title 'P_{21}' with lines,  \
	    
	       

