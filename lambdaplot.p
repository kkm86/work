# Gnuplot script file for plotting data in file "infty~.dat"
      # This file is called   lambdaplot.p
  
   set term aqua enhanced font "Times, 18" dashed   # turn on enhanced text mode
   
   
#set termoption enhanced
   set encoding utf8
   set minussign
   set   autoscale                        # scale axes automatically
   unset log                              # remove any log-scaling
   unset label                            # remove any previous labels
      unset arrow
      unset xtic
      unset title
      unset zeroaxis
      
   set xzeroaxis
      set key right top
      set style line 1 dt 1 lc rgb 'red' pt 13 
      set style line 2 dt 5 lc rgb 'blue' pt 2 
      set style line 3 dt 3 lc rgb '#008B8B' pt 8
      set style line 4 dt 4 lc rgb '#696969' pt 10
      set style line 5 dt 1 lc rgb '#00CED1' pt 11
#set arrow from 1.0,-1.0125145 to 20000.0,-1.0125 nohead lw 1 dt 2 lc 'black'
      set arrow from 55.0,-0.2 to 55.0,0.2 nohead lw 1 dt 2 lc 'black'
      #set title font "TMR,20"
      set xlabel font "Times,20"
      set ylabel font "Times,20"         
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically
      set key font "Times,18"
      set xlabel "{/Symbol r}/{/Times-Italic r}_0"
      set ylabel "{/Symbol x}({/Symbol r})"
#set ylabel "2{/Symbol mr}² {/Times-Italic U}_{/Symbol n} ({/Symbol r}) + 1/4"
      set log x
      set xr [1.0:20000.0]
      set yr [-0.2:0.2]
      plot "inftyneg.dat" using 2:((10**8)*(($3)-0.25)/(2*92228.8*($2**2))) title "{/Times-Italic=18 a}_1" with lines ls 1,  \
	       "inftypos.dat" using 2:((10**8)*(($3)-0.25)/(2*92228.8*($2**2))) title"{/Times-Italic=18 a}_2" with lines ls 2,  \
		  "inftypos.dat" using 2:((10**8)*(-1.01252-0.25)/(2*92228.8*($2**2))) title"{/Times-Italic=18 a}_2" with lines ls 2,  \
     
		 

			      
