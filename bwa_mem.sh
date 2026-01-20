CPU=8
REF=Xiphophorus_maculatus.X_maculatus-5.0-male.dna_sm.toplevel.fa
bwa=/usr/bin/bwa

#step 1  index

bwa index $REF

#step 2 run bwa mem & samtools partly

for sR1 in *.paired_1.fq.gz;do
        sBase=`dirname $sR1`
        sStem=`basename $sR1`
        sName=${sStem/.paired_1.fq.gz/}
        sR2=$sBase/$sName.paired_2.fq.gz

        ( $bwa mem -t $CPU $REF $sR1 $sR2 \
                |samtools view --threads 2 -F 0x4 -u -|samtools sort -m 10g --threads 8 -o $sName.bam ) > $sName.log 2>&1 &

done
wait
