# Gnuplot script file for plotting data in file "result.dat"
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
      set ytic auto
      set ztic auto
                                             # set ytics automatically
      set key font ",12"
      set title "U(R)"
	 set xlabel "R/a"
      set ylabel
	 set log x
	 set xr [0.01:50.00]
	 set yr [-5000.0:4.0]
#set yr [-2.0:2.37]
	 plot "result2.dat" using 2:($3*2) title '3/4'  with linespoint,  \
	    "result3.dat" using 2:($3*2) title '-1/4'  with linespoint,  \
     "result4.dat" using 2:($3*2) title '+15/4'  with linespoint,  \
	       "result5.dat" using 2:($3*2) title 'ny'  with linespoint,  \

