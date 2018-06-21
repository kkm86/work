# Gnuplot script file for plotting data in file "threebodypot.dat"
      # This file is called   wave.p
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
	 set ylabel "U(R)+V_{trap}(R) (osc.units)"
	 set log x
	 set xr [0.01:2.0]
	 set yr [-5000.0:500.0]
	 plot "threebodypot.dat" using 2:3 title 'v=0'  with lines,  \
	    "threebodypot.dat" using 2:4 title 'v=1'  with lines,  \
	       "threebodypot.dat" using 2:5 title 'v=2'  with lines,  \
		  "threebodypot.dat" using 2:6 title 'v=3'  with lines,  \
		     "wave.dat" using 2:3 title 'N²'  with lines,  \
#	"wave.dat" using ($2/731):3 title 'Squarred norm'  with lines,  \
	  
	       

