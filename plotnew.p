# Gnuplot script file for plotting data in file "effectivepot0_4.dat"
      # This file is called   plotnew.p
set termopt enhanced    # turn on enhanced text mode
set encoding utf8
      
      set   autoscale                        # scale axes automatically
      unset log                              # remove any log-scaling
      unset label                            # remove any previous labels
      unset arrow
      unset xtic
      unset title
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
#set title ""
         set xlabel "{/Symbol r}"
         set ylabel "2{/Symbol mr}² W_{/Symbol n} ({/Symbol r}) + 1/4"
         set log x
         set xr [1.0:1000000.0]
         set yr [-5.0:5.0]
         plot "effectivepot0_45_N5.dat" using 2:3 title "{/Times:Italic=12 a{/TMR:Normal=10 =1018}}"  with points ls 3,  \
	    "effectivepot0_45_N20.dat" using 2:3 notitle with points ls 3,  \
	       "effectivepot0_42_p0.dat" using 2:3 title "{/Times:Italic=12 a{/TMR:Normal=10 =5654}}" with points ls 1,  \
		  "effectivepot0_42_p1.dat" using 2:3 notitle with points ls 1,  \
		     "effectivepot0_42.dat" using 2:3 notitle with points ls 1,  \
			"effectivepot0_414_p0.dat" using 2:3 title "{/Times:Italic=12 a{/TMR:Normal=10 =291011}}"  with points ls 2,  \
			   "effectivepot0_414_p0.dat" using 2:4 notitle with points ls 2,  \
			      "effectivepot0_414_p1.dat" using 2:3 notitle with points ls 2,  \
				 "effectivepot0_414_N5.dat" using 2:3 notitle with points ls 2,  \
