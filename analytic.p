# Gnuplot script file for plotting data in file "lami~.dat"
      # This file is called   analytic.p
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
#set log x
	 set xr [1.0:60.0]
         set yr [0.0:4.5]
         plot "lami20nyy.dat" using 2:3 title 'N_{/Symbol q}=20'  with lines,  \
	    "lami30nyy.dat" using 2:3 title 'N_{/Symbol q}=30'  with lines,  \
	       "lami40nyy.dat" using 2:3 title 'N_{/Symbol q}=40'  with lines,  \
	       "lamii.dat" using 1:2 title 'Analytic'  with lines,  \

	    
