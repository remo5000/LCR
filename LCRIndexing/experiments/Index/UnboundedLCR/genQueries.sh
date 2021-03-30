#!/bin/sh

# TODO set pattern here
datasets=`ls experiments/graphs/real/stringFc.edge | sed 's:\.edge::g'`;
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

    tmp=`mktemp`
    if [ ! -f "${ds}.querylog" ];
    then
        echo "L=${L}";
        echo "build/default/genQuery ${ds}.edge 3 1000 ${L1} ${L2} ${L3} > ${tmp}";
        build/default/genQuery ${ds}.edge 3 1000 ${L1} ${L2} ${L3} > ${tmp}
        echo "cp ${tmp} ${ds}.querylog";
        cp ${tmp} ${ds}.querylog
        echo "rm ${tmp}";
        rm ${tmp}
    fi
done;
