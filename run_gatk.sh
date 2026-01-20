#!/bin/bash
#do gatk.sh for each sample.bam file

for sSample in *.bam; do
        sBase=`basename $sSample`
        sName=${sBase/.bam/}

        source gatk.sh $sSample &
done
wait
