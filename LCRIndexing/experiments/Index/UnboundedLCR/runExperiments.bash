#!/bin/bash
#datasets=`find ./experiments/graphs/ | egrep '\.edge$' | sed 's:\.edge::g' | egrep -i '(V5kD(2|3|4|5)L8exp)'`;
# datasets=`find ./experiments/graphs/ | egrep '\.edge$' | sed 's:\.edge::g' | egrep 'L8exp' | egrep 'V5kD'`;
datasets=`find ./experiments/graphs/ | egrep '\.edge$' | sed 's:\.edge::g'`;
numberOfMethods=6
numberOfQuerySets=1

#datasets=`find ./experiments/graphs/ | egrep '\.edge$' | sed 's:\.edge::g' | egrep -i '(notre|bio)'`;


for ds in $datasets;
do
        dsp=`echo $ds | sed 's/.*\///'`;
        echo "build/default/runExperiment "$ds".edge $numberOfQuerySets "$1"results-"$dsp".csv > "$1"output-"$dsp".txt";
        build/default/runExperiment $ds".edge" $numberOfQuerySets $1"results-"$dsp".csv" > $1"output-"$dsp".txt";
done;
