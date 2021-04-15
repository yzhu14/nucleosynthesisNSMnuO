set term postscript enhanced color font 'Helvetica,24'
#set output '| epstopdf --filter --outfile=nu_Surf567eff0_6005_2069.pdf'
#set output 'AMfluxFrensel2a0.35abA.eps'
set output 'nonuAB.eps'
set yrange [1e-12:5e-2]
set xrange [0:140]
#set ylabel "P(s)" font 'Helvetica,24'
#set xlabel "s" font 'Helvetica,24'
set border lw 5
set logscale y
set format y "10^{%L}"
#set logscale x
#set yrange [0.00000001:0.001]
#set key width -3 spacing 1. height .1
set key top right horizontal font "Helvetica-Bold,12"
set style line 100 lt 1 lc rgb "black" lw 5
set style line 1 lt 1 lc rgb "#a6cee3" lw 5
set style line 2 lt 1 lc rgb "#1f78b4" lw 5
set style line 3 lt 1 lc rgb "#b2df8a" lw 5
set style line 4 lt 1 lc rgb "#33a02c" lw 2
set style line 5 lt 1 lc rgb "#fb9a99" lw 2
set style line 6 lt 1 lc rgb "#e31a1c" lw 2
set style line 7 lt 1 lc rgb "#fdbf6f" lw 2
set style line 8 lt 1 lc rgb "#ff7f00" lw 2
set style line 9 lt 1 lc rgb "#cab2d6" lw 2
set style line 10 lt 1 lc rgb "#6a3d9a" lw 2
set style line 11 lt 3 lc rgb "#a6cee3" lw 5
set style line 12 lt 3 lc rgb "#1f78b4" lw 5
set style line 13 lt 3 lc rgb "#b2df8a" lw 5
set style line 14 lt 3 lc rgb "#33a02c" lw 2
set style line 15 lt 3 lc rgb "#fb9a99" lw 2
set style line 16 lt 3 lc rgb "#e31a1c" lw 2
set style line 17 lt 3 lc rgb "#fdbf6f" lw 2
set style line 18 lt 1 lc rgb "#ff7f00" lw 2
set style line 19 lt 3 lc rgb "#cab2d6" lw 2
set style line 20 lt 3 lc rgb "#6a3d9a" lw 2
set style line 21 lt 1 lc rgb "red" lw 5
set style line 22 lt 1 lc rgb "blue" lw 5
set style line 23 lt 1 lc rgb "purple" lw 5
set style line 24 lt 1 lc rgb "orange" lw 5
set style line 25 lt 3 lc rgb "cyan" lw 5
set style line 26 lt 3 lc rgb "orange" lw 2
set style line 27 lt 3 lc rgb "blue" lw 2
set style line 28 lt 1 lc rgb "green" lw 5
set style line 29 lt 1 lc rgb "yellow" lw 2
set style line 30 lt 1 lc rgb "black" lw 3
set style line 31 lt 3 lc rgb "red" lw 5


plot "/Users/yzhu14/Research/nucleosynthesisNSMnuO/albino/oldtests/testTracerofDirk/finab/finab_57221.dat" using 1:($5/100) with point lw 5 ti "57221",\
"nonu57221abA" using 1:2 with points lw 5 ti "nonu57221",\
#"noncap57221abA" using 1:2 with points lw 5 ti "noncap57221",\

plot "/Users/yzhu14/Research/nucleosynthesisNSMnuO/albino/oldtests/testTracerofDirk/finab/finab_66394.dat" using 1:($5/100) with point lw 5 ti "66394",\
"nonu66394abA" using 1:2 with points lw 5 ti "nonu66394",\
#"noncap66394abA" using 1:2 with points lw 5 ti "noncap66394"

plot "/Users/yzhu14/Research/nucleosynthesisNSMnuO/albino/oldtests/testTracerofDirk/finab/finab_78323.dat" using 1:($5/100) with point lw 5 ti "78323",\
"nonu78323abA" using 1:2 with points lw 5 ti "nonu78323",\
#"noncap78323abA" using 1:2 with points lw 5 ti "noncap78323",\

plot "/Users/yzhu14/Research/nucleosynthesisNSMnuO/albino/oldtests/testTracerofDirk/finab/finab_79049.dat" using 1:($5/100) with point lw 5 ti "79049",\
"nonu79049abA" using 1:2 with points lw 5 ti "nonu79049",\
#"noncap79049abA" using 1:2 with points lw 5 ti "noncap79049",\

plot "/Users/yzhu14/Research/nucleosynthesisNSMnuO/albino/oldtests/testTracerofDirk/finab/finab_80224.dat" using 1:($5/100) with point lw 5 ti "80224",\
"nonu80224abA" using 1:2 with points lw 5 ti "nonu80224",\
#"noncap80224abA" using 1:2 with points lw 5 ti "noncap80224",\

plot  "/Users/yzhu14/Research/nucleosynthesisNSMnuO/albino/oldtests/testTracerofDirk/finab/finab_61778.dat" using 1:($5/100) with point lw 5 ti "61778",\
"nonu61778abA" using 1:2 with points lw 5 ti "nonu61778",\
#"noncap61778abA" using 1:2 with points lw 5 ti "noncap61778"


