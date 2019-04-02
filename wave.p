# Gnuplot script file for plotting data in file "effectivepot.dat"
      # This file is called   wave.p
 set termopt enhanced    # turn on enhanced text mode
   
      
   
      set   autoscale                        # scale axes automatically
      unset log                              # remove any log-scaling
      unset label                            # remove any previous labels
      unset arrow
      unset xtic
      unset title
      set encoding utf8
      set xzeroaxis
#set title font "Times,20"
      set xlabel font "Times,15"
      set ylabel font "Times,15"         
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically
      set key font ",12"
#set title "U({/Symbol r})"
         set xlabel "{/Symbol r}/{/Helvetica a}"
         set ylabel "U({/Symbol r}) [10^{-8} a.u.]"
         set log x
         set xr [0.03:100.4]
         set yr [-6.0:20.0]
         plot "effective.dat" using 2:3 title '{/Symbol n}=0'  with lines,  \
	    "effective.dat" using 2:4 title '{/Symbol n}=4'  with lines,  \
	       "effective.dat" using 2:5 title '{/Symbol n}=6'  with lines,  \
		  "effective.dat" using 2:6 title '{/Symbol n}=0 gegen'  with lines,  \
		     "effective.dat" using 2:7 title '{/Symbol n}=4  gegen'  with lines,  \
			"effective.dat" using 2:8 title '{/Symbol n}=6 gegen'  with lines,  \
	       

	  
	       

