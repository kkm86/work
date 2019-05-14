# Gnuplot script file for plotting data in file "infty~.dat"
      # This file is called   lambdaplot.p
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
      set style line 1 dt 1 lc rgb 'red' pt 13 
      set style line 2 dt 5 lc rgb 'blue' pt 2 
      set style line 3 dt 3 lc rgb '#008B8B' pt 8
      set style line 4 dt 4 lc rgb '#696969' pt 10
      set style line 5 dt 1 lc rgb '#00CED1' pt 11
      set arrow from 1.0,-1.0125145 to 10000.0,-1.0125 nohead lw 1 dt 2 lc 'black'
      #set title font "TMR,20"
      set xlabel font "Calibri,20"
      set ylabel font "Calibri,20"         
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically
      set key font "Calibri,16"
      set xlabel "{/Symbol r} [a.u.]"
      set ylabel "2{/Symbol mr}² U_{/Symbol n} ({/Symbol r}) + 1/4"
      set log x
      set xr [1.0:20000.0]
      set yr [-5.0:4.1]
      plot"inftyneg.dat" using 2:3 title "{/Times-Italic=18 a_1" with lines ls 1,  \
	       "inftypos.dat" using 2:3 title"{/Times-Italic=18 a_2" with lines ls 2,  \
#"lambda3p_tot.dat" using 2:3 title"{/Times:Italic=12 a{/TMR:Normal=10 =30632}}" with lines ls 2,  \
			"lambda4p_tot.dat" using 2:3 title "{/Times:Italic=12 a{/TMR:Normal=10 =291010}}" with lines ls 1,  \
			   "lambda4p_tot.dat" using 2:(-1/(87*1836.15*(1018**2))) title "{/Times:Italic=12 E{/TMR:Normal=10 =}}" with lines ls 1,  \
