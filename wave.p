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
         set xr [0.0003:100.4]
         set yr [-16.0:5.0]
         plot "effectivepot.dat" using 2:3 title '{/Symbol n}=0, neg a'  with lines,  \
	    "effectivepotpos.dat" using 2:3 title '{/Symbol n}=0 pos a'  with lines,  \
	       "effectivepotneg.dat" using 2:3 title '{/Symbol n}=0 neg large'  with lines,  \
		  "effectivepotneg2.dat" using 2:3 title '{/Symbol n}=0 neg largest'  with lines,  \
	       

	  
	       

