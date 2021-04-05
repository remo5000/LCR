#!/bin/sh

set -x 

# TODO set pattern here
datasets=`ls experiments/graphs/real/web-Stanford.edge | sed 's:\.edge::g'`;
for ds in $datasets;
do
    dsp=`echo $ds | sed 's/.*\///'`;
    L=`cat ${ds}.edge | awk '{ if( $3 in used ){ } else { used[$3]=1; print $3 } }' | wc -l`;
    # # TODO set number of labels here
    # L1=1;
    # L2=2;
    # L3=2;
    L1=`expr $L / 4`;
    L2=`expr $L / 2`;
    L3=`expr $L - 2`;
    NUM_QUERIES=10;
	

    tmp=`mktemp`
    if [ ! -f "${ds}.querylog" ];
    then
        build/default/genQuery ${ds}.edge 3 ${NUM_QUERIES} ${L1} ${L2} ${L3} > ${tmp};
        cp ${tmp} ${ds}.querylog;
        rm ${tmp};
    fi
done;
