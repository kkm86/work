# Gnuplot script file for plotting data in file "N*tot.dat"
      # This file is called   universal.p
   set term aqua enhanced font "Calibri, 20" dashed   # turn on enhanced text mode
   
      
   
   set   autoscale                        # scale axes automatically
   unset log                              # remove any log-scaling
   unset label                            # remove any previous labels
      unset arrow
      unset xtic
      unset title
      unset zeroaxis
      set encoding utf8
#set xzeroaxis
      set key left bottom
      set style line 1 dt 1 lc rgb '#C71585' pt 13 
      set style line 2 dt 5 lc rgb '#000080' pt 9 
      set style line 3 dt 3 lc rgb '#008B8B' pt 8
      set style line 4 dt 4 lc rgb '#696969' pt 10
      set style line 5 dt 1 lc rgb '#00CED1' pt 11
      set arrow from 1.0,-1.0125 to 20000.0,-1.0125 nohead lw 1 dt 2 lc 'black'
      #set title font "TMR,20"
      set xlabel font "Calibri,20"
      set ylabel font "Calibri,20"         
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically
      set key font "Calibri,16"
      set xlabel "{/Symbol r} [a.u.]"
      set ylabel "2{/Symbol mr}² {/:Italic=18 W}_{/Symbol n} ({/Symbol r}) + 1/4"
      set log x
      set xr [1.0:20000.0]
      set yr [-5.0:4.1]
      plot"N1tot.dat" using 2:3 title "{/Times-Italic=16 a {/Calibri:Normal=16 = -2386}} a.u." with lines ls 4,  \
	 "N2tot.dat" using 2:3 title"{/Times-Italic=16 a {/Calibri:Normal=16 = -8721}} a.u." with lines ls 3,  \
	    "N3tot.dat" using 2:3 title"{/Times-Italic=16 a {/Calibri:Normal=16 = -38612}} a.u." with lines ls 2,  \
	       "N4tot.dat" using 2:3 title "{/Times-Italic=16 a {/Calibri:Normal=16 = -2702020}} a.u." with lines ls 1,  \
			   
