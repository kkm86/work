# Gnuplot script file for plotting data in file "U~.dat"
      # This file is called   threebodypot.p
   set term aqua enhanced font "Calibri, 18" dashed   # turn on enhanced text mode
 
   
      
   
      set   autoscale                        # scale axes automatically
      unset log                              # remove any log-scaling
      unset label                            # remove any previous labels
      unset arrow
      unset xtic
      unset title
      unset zeroaxis 
      set encoding utf8
      set xzeroaxis linetype 1 dashtype 3 lc 'black'
#set title font "Times,20"
      set xlabel font "Calibri,20"
      set ylabel font "Calibri,20"         
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically
      set key font "Calibri,16"
      set key right bottom
#set title "U({/Symbol r})"
         set xlabel "{/Symbol r}/|{/Times-Italic=20 a}|"
         set ylabel "W_{{/Symbol n}}({/Symbol r}) [10^{-10} a.u.]"
	 set log x
	 set xr [0.044:20.0]
         set yr [-50.0:15.0]
         plot "Upos.dat" using 2:3 title '{/Symbol n}=0'  with lines,  \
	    "Upos.dat" using 2:4 title '{/Symbol n}=1' with lines,  \
	       "Upos.dat" using 2:5 title '{/Symbol n}=2'  with lines,  \
		  "Upos.dat" using 2:6 title '{/Symbol n}=3'  with lines,  \
		     "Upos.dat" using 2:7 title '{/Symbol n}=4'  with lines,  \
			"Upos.dat" using 2:8 title 'E_{2b}'  with lines dt 2,  \
	    
