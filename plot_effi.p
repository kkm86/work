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
      set title "Two-body Potential"
      set xlabel "Bohr radius"
      set ylabel 
	 set xr [0.0:10.0]
	 set yr [-100.6:0.5]
	 plot    "result.dat" using 2:3 title 'TBT'  with linespoint,  \
	  
	       

