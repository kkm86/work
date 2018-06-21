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
      set title "U(R)"
	 set xlabel "R(osc.units)"
	 set ylabel "Squared norm"
	 set log x
	 set xr [0.01:2.0]
	 set yr [0.9:1.3]
	 plot "wave.dat" using 2:3 notitle with lines,  \
	    
	       

