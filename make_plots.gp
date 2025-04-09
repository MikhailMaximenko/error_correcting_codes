unset title
set xlabel 'E_b / N_0, dB'
set ylabel 'error rate'
set logscale y
set ytics format '10^{%T}'
set xtics 0.5
set mytics 10
set grid xtics mytics
set terminal png size 1024,768 enhanced
set output 'results.png'
plot for [i=1:|ARGV|] ARGV[i] using ($5):(1 - $4) with linespoints title ARGV[i] noenhanced