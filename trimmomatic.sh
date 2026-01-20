trimjar=/data/software/Trimmomatic-0.39/trimmomatic-0.39.jar
adapterfa=/data/software/Trimmomatic-0.39/adapters/bgi_adapter.fa

for sR1 in *_1.fq.gz;do
        sStem=`basename $sR1`
        sStem=${sStem/%_1.fq.gz}
        sPath=`dirname $sR1`
        sR2=$sPath/${sStem}_2.fq.gz
        java -jar $trimjar PE -threads 16 -phred33 $sR1 $sR2 $sStem.paired_1.fq.gz $sStem.unpaired_1.fq.gz  $sStem.paired_2.fq.gz $sStem.unpaired_2.fq.gz  ILLUMINACLIP:$adapterfa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > trim.$sStem.log 2>&1 &
done
wait
