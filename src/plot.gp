reset
set auto
set xlabel "alpha"
set ylabel "|| x - x^*||_1 / N"
set arrow from alphac,0 to alphac,0.9 nohead
p "<awk '$1==100' results.dat" u 2:3:4 title "N=100" w e,\
"<awk '$1==200' results.dat" u 2:3:4 title "N=200" w e,\
"<awk '$1==500' results.dat" u 2:3:4 title "N=500" w e

alphac=.6455722919969774
