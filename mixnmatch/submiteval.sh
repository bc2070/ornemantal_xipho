here=`pwd`
for i in sims/*/;do
	cd $here;
	cd $i;
	sbatch evaluate.sh
done
