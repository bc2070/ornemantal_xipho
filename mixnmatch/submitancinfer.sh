here=`pwd`
for i in sims/*/;do
	cd $here;
	cd $i;
	php $here/thinaims.php simulated_AIMs_for_AncestryHMM simulated_parental_counts_for_AncestryHMM simulated_AIMs_for_AncestryHMM_thinned simulated_parental_counts_for_AncestryHMM_thinned
	sbatch ancinfer.sbatch
done

