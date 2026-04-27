for i in target_region2/*;do
	sBase=`basename $i`
        sName=${sBase/%.bed}
	python3 02.match_id.py Xiphophorus_maculatus.longest.isoform.mRNA.gff3 $i ${sName}_matched.bed
	python3 03.cal_d.py ${sName}_matched.bed dN_stats_genes.txt Xiphophorus_maculatus.longest.isoform.mRNA.gff3 ${sName}_cal.bed
	awk -F'\t' '$4 != "."' ${sName}_cal.bed > ${sName}_cal_filtered.bed
	rm ${sName}_matched.bed ${sName}_cal.bed
done
wait
