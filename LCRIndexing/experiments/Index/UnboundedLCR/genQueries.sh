#!/bin/sh

set -x 

datasets=`find experiments/graphs/ | egrep '\.edge$' | sed 's:\.edge::g'`;
for ds in $datasets;
do
    dsp=`echo $ds | sed 's/.*\///'`;
    L=`cat ${ds}.edge | awk '{ if( $3 in used ){ } else { used[$3]=1; print $3 } }' | wc -l`;
    L1=`expr $L / 4`;
    L2=`expr $L / 2`;
    L3=`expr $L - 2`;
    NUM_QUERIES=1000;
	
    tmp=`mktemp`
    if [ ! -f "${ds}.querylog" ];
    then
        build/default/genQuery ${ds}.edge 3 ${NUM_QUERIES} ${L1} ${L2} ${L3} > ${tmp};
        cp ${tmp} ${ds}.querylog;
        rm ${tmp};
    fi
done;
