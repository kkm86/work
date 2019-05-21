# Gnuplot script file for plotting data in file "lami~.dat"
      # This file is called   faddeev.p
   set term aqua enhanced font "Times, 18" dashed   # turn on enhanced text mode
 
   
      
      set minussign
      set   autoscale                        # scale axes automatically
      unset log                              # remove any log-scaling
      unset label                            # remove any previous labels
      unset arrow
      unset xtic
      unset title
      set encoding utf8
      set xzeroaxis
      set arrow from -4.0,-1.0125 to 4.0,-1.0125 nohead lw 1 dt 2 lc 'black'
#set title font "Times,20"
      set xlabel font "Times,20"
      set ylabel font "Times,20"        
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically
      set key font "Times,18"
      set key right bottom
#set title "U({/Symbol r})"
         set xlabel "log_{10}({/Symbol r}/ |{/Times-Italic=20 a}|)"
         set ylabel "{/Symbol n}_n({/Symbol r}/ |{/Times-Italic=20 a}|)"
#set log x
	 set xr [-4.0:4.0]
         set yr [-15.0:20.0]
         plot "lami1.dat" using 1:2 title '{/Times-Italic=20 a} > 0'  with lines dt 1,  \
	    "lami2.dat" using 1:2 title '{/Times-Italic=20 a} < 0'  with lines dt 1,  \
	       "lami3.dat" using 1:2 title '{/Times-Italic=20 a} > 0'  with lines dt 4,  \
	       

	    
