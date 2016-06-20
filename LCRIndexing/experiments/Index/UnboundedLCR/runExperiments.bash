#!/bin/bash
#datasets=`find ./experiments/graphs/ | egrep '\.edge$' | sed 's:\.edge::g' | egrep -i '(/V25kD5L8exp)' | sort -d`; # N/50, b=20
#datasets=`find ./experiments/graphs/ | egrep '\.edge$' | sed 's:\.edge::g' | egrep -i '(jmd|citeSeer|NotreDame|soc-sign|Slashdot|V125|Twitter)'`; # N/100, b=20
#datasets=`find ./experiments/graphs/ | egrep '\.edge$' | sed 's:\.edge::g' | egrep -i '(V625k|usPatents|webBerk|webGoogle)'`; # N/500, (b=15,b=30)
datasets=`find ./experiments/graphs/ | egrep '\.edge$' | sed 's:\.edge::g' | egrep -i '(webStanford)'`;

for ds in $datasets;
do
    if [[ "`cat $ds".edge" | wc -l`" -le $2 &&  "`cat $ds".edge" | wc -l`" -gt $3 ]]; then
        dsp=`echo $ds | sed 's/.*\///'`;
        echo "build/default/runExperiment "$ds".edge 3 "$1"results-"$dsp".csv > "$1"output-"$dsp".txt";
        build/default/runExperiment $ds".edge" 3 $1"results-"$dsp".csv" > $1"output-"$dsp".txt";
    fi
done;
