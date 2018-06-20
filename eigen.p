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
	 set xlabel "R[Osc. units]"
	 set ylabel "V_{trap}+U_{1}[osc. units]"
	 set log x
	 set xr [1:1000.0]
	 set yr [-6.0:2.0]
#set yr [-2.0:2.37]
	 plot "result10.dat" using 2:3 notitl with lines,  \
	    "result10.dat" using 2:4 notitle  with lines,  \
	       "result10.dat" using 2:5 title '2-b energy' with lines,  \
#  "result10.dat" using 2:5 title '2-b energy' with lines,  \
		
