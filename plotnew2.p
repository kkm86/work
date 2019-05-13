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
      set key right bottom
      set style line 1 lc rgb '#C71585' pt 13 
      set style line 2 lc rgb '#000080' pt 9 
      set style line 3 lc rgb '#008B8B' pt 8
      set arrow from 1.0,-1.0125 to 1000000.0,-1.0125 nohead lw 1 dt 2 lc 'grey'
      #set title font "TMR,20"
      set xlabel font "TMR,15"
      set ylabel font "TMR,15"         
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically
      set key font ",12"
      set log x
      set xr [1.0:1000000.0]
      set yr [-5.0:5.0]
      plot "effectivepot0_4.dat" using 2:3 notitle with points ls 1,  \
	 "effectivepot0_4_p0.dat" using 2:3 notitle  with points ls 1,  \
	    "effectivepot0_4_p3.dat" using 2:3 notitle  with points ls 1,  \
	       "effectivepot0_4_p.dat" using 2:3 notitle  with points ls 1,  \
		  "effectivepot0_4_p2.dat" using 2:3 title "{/Times:Italic=12 a{/TMR:Normal=10 =-2386}}"  with points ls 1,  \
		     "effectivepot0_41.dat" using 2:3 notitle  with points ls 2,  \
			"effectivepot0_41_p0.dat" using 2:3 title "{/Times:Italic=12 a{/TMR:Normal=10 =-8721}}"  with points ls 2,  \
			   "effectivepot0_413_N10.dat" using 2:3 notitle  with points ls 3,  \
			      "effectivepot0_413_N20.dat" using 2:3 notitle  with points ls 3,  \
				 "effectivepot0_413_N30.dat" using 2:3 notitle   with points ls 3,  \
				    "effectivepot0_413_N60.dat" using 2:3 title "{/Times:Italic=12 a{/TMR:Normal=10 =-38612}}"  with points ls 3,  \
				       "effectivepot0_413_N80.dat" using 2:3 title "{/Times:Italic=12 a{/TMR:Normal=10 =-38612}}"  with points,  \
					  
   
			   
