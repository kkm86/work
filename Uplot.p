# Gnuplot script file for plotting data in file "U~.dat"
      # This file is called   threebodypot.p
   set term aqua enhanced font "Times, 18" dashed   # turn on enhanced text mode
 
   
      
   
      set   autoscale                        # scale axes automatically
      unset log                              # remove any log-scaling
      unset label                            # remove any previous labels
      unset arrow
      unset xtic
      unset title
      unset zeroaxis 
      set encoding utf8
      set xzeroaxis linetype 1 dashtype 3 lc 'black'
	 set arrow from 0.218546,-50. to 0.218546,15.0 nohead lw 1 dt 2 lc 'green'
	 set arrow from 1.0,-50. to 1.0,15.0 nohead lw 1 dt 2 lc 'green'
#set title font "Times,20"
      set xlabel font "Times,20"
      set ylabel font "Times,20"         
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically
      set key font "Times,18"
      set key right bottom
#set title "U({/Symbol r})"
         set xlabel "{/Symbol r}/|{/Times-Italic=20 a}|"
         set ylabel "{/Times-Italic=20 W}_{{/Symbol n}}({/Symbol r}) [10^{-10} a.u.]"
	 set log x
	 set xr [0.044:20.0]
         set yr [-50.0:15.0]
         plot "Upos.dat" using 2:3 title '{/Symbol n} = 0'  with lines,  \
	    "Upos.dat" using 2:4 title '{/Symbol n} = 1' with lines,  \
	       "Upos.dat" using 2:5 title '{/Symbol n} = 2'  with lines,  \
		  "Upos.dat" using 2:6 title '{/Symbol n} = 3'  with lines,  \
		     "Upos.dat" using 2:7 title '{/Symbol n} = 4'  with lines,  \
			"Upos.dat" using 2:8 title '{/Times-Italic E}_{2b}'  with lines dt 2,  \
			   
