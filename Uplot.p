# Gnuplot script file for plotting data in file "U~.dat"
      # This file is called   threebodypot.p
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
      set key right bottom
#set title "U({/Symbol r})"
         set xlabel "{/Symbol r}/{/Times:Italic=15 |a|}"
         set ylabel "W_{{/Symbol n}}({/Symbol r}) [10^{-10} a.u.]"
	 set log x
	 set xr [0.044:20.0]
         set yr [-6.0:4.0]
         plot "Uneg.dat" using 2:3 title '{/Symbol n}=0'  with lines,  \
	    "Uneg.dat" using 2:4 title '{/Symbol n}=1' with lines,  \
	       "Uneg.dat" using 2:5 title '{/Symbol n}=2'  with lines,  \
		  "Uneg.dat" using 2:6 title '{/Symbol n}=3'  with lines,  \
		     "Uneg.dat" using 2:7 title '{/Symbol n}=4'  with lines,  \
#	"Uneg.dat" using 2:8 title 'E_{2b}'  with lines dt 2,  \
	    
