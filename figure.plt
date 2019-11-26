set xrange [0:70]
set yrange [0:70]
set term png
set output "result.png"

p 'out.txt' using 2:3 with linespoints notitle, 'out.txt' using 2:3:(sprintf("[%d]", $1)) with labels offset .5,.5 font ",7" notitle