# Gnuplot script file for plotting data in file "effectivepot0_4.dat"
      # This file is called   plotnew2.p
 set termopt enhanced    # turn on enhanced text mode
   
      
   
      set   autoscale                        # scale axes automatically
      unset log                              # remove any log-scaling
      unset label                            # remove any previous labels
      unset arrow
      unset xtic
      unset title
      set encoding utf8
      set xzeroaxis
	 set style line 1 lc rgb 'blue' pt 5 #square
	 set style line 2 lc rgb 'green' pt 7 #circle
	 set style line 3 lc rgb 'red' pt 9 #triangle
	 set arrow from 1.0,-1.0125 to 1000000.0,-1.0125 nohead
#set title font "Times,20"
      set xlabel font "Times,15"
      set ylabel font "Times,15"         
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically
      set key font ",12"
#set title "U({/Symbol r})"
         set xlabel "{/Symbol r}"
	 set ylabel "2{/Symbol m}{/Symbol r}²*W({/Symbol r}) + 1/4"
         set log x
         set xr [1.0:1000000.0]
         set yr [-3.0:5.0]
         plot "effectivepot0_4.dat" using 2:3 notitle with points ls 1,  \
		     "effectivepot0_4_p0.dat" using 2:3 notitle  with points ls 1,  \
			"effectivepot0_4_p3.dat" using 2:3 notitle  with points ls 1,  \
			   "effectivepot0_4_p.dat" using 2:3 notitle  with points ls 1,  \
			   "effectivepot0_4_p2.dat" using 2:3 title 'a=-2385.88'  with points ls 1,  \
				    "effectivepot0_41.dat" using 2:3 notitle  with points ls 2,  \
				       "effectivepot0_41_p0.dat" using 2:3 title 'a=-8720.93'  with points ls 2,  \
					  "effectivepot0_413_N20.dat" using 2:3 notitle  with points ls 3,  \
					  "effectivepot0_413_N5.dat" using 2:3 notitle  with points ls 3,  \
					     "effectivepot0_413_N20long.dat" using 2:3 title 'a=-38612.2'  with points ls 3,  \
					  
   
			   
