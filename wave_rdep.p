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
	 set xr [0.01:50.00]
	 set yr [-5000.0:500.0]
#set yr [-2.0:2.37]
	 plot "result5.dat" using 2:3 title 'Blume'  with lines,  \
	    "result6.dat" using 2:3 notitle  with lines,  \

