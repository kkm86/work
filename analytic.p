# Gnuplot script file for plotting data in file "lami~.dat"
      # This file is called   analytic.p
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
#set title font "Times,20"
      set xlabel font "Times,20"
      set ylabel font "Times,20"        
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically
      set key font "Times,18"
      set key right bottom
#set title "U({/Symbol r})"
         set xlabel "{/Symbol r}/ |{/Times-Italic=20 a}|"
         set ylabel "2{/Symbol mr}² {/Times-Italic=20 W}_{/Symbol n} ({/Symbol r}) + 1/4"
#set log x
	 set xr [1.0:60.0]
         set yr [0.0:4.5]
         plot "lami20nyy.dat" using 2:3 title 'N_{/Symbol q} = 20'  with lines dt 4,  \
	    "lami30nyy.dat" using 2:3 title 'N_{/Symbol q} = 30'  with lines dt 5,  \
	       "lami40nyy.dat" using 2:3 title 'N_{/Symbol q} = 40'  with lines,  \
	       "lamii.dat" using 1:2 title 'Analytic'  with lines lc 'black',  \

	    
