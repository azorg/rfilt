#!/usr/bin/env gnuplot

set grid
unset border
plot "data.txt" using ($1):($2) title "in"     with lines, \
     "data.txt" using ($1):($3) title "out"    with lines, \
     "data.txt" using ($1):($4) title "v"      with lines, \
     "data.txt" using ($1):($5) title "a"      with lines, \
     "data.txt" using ($1):($6) title "r"      with lines, \
     "data.txt" using ($1):($7) title "x_stop" with lines
pause -1
     

