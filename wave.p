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
         set xr [1.0:1000000.0]
         set yr [-8.0:5.0]
         plot "effectivepot0_45.dat" using 2:3 title '{/Symbol n}=0'  with lines,  \
	    "effectivepot0_45_2.dat" using 2:3 title '{/Symbol n}=1'  with lines,  \
	       "effectivepot0_65.dat" using 2:3 title '{/Symbol n}=0'  with lines,  \
		  "effectivepot0_45.dat" using 2:6 title '{/Symbol n}=1'  with lines,  \
		     "effectivepot0_4.dat" using 2:3 title '{/Symbol n}=0'  with lines,  \
		     "effectivepot0_4_2.dat" using 2:3 title '{/Symbol n}=1'  with lines,  \
			"effectivepot0_65_2.dat" using 2:3 title '{/Symbol n}=0'  with lines,  \
			"effectivepot0_65_p.dat" using 2:3 title '{/Symbol n}=1'  with points,  \
			   "effectivepot0_65_p2.dat" using 2:3 title '{/Symbol n}=1'  with points,  \
			      "effectivepot0_4_p3.dat" using 2:3 title '{/Symbol n}=1'  with points,  \
				 "effectivepot0_4_p.dat" using 2:3 title '{/Symbol n}=1'  with points,  \
				    "effectivepot0_4_p0.dat" using 2:3 title '{/Symbol n}=1'  with points,  \
		
	       

	  
	       

