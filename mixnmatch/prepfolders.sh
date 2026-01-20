arrDrift=( 0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 )
sOutDir=sims

for nDrift in "${arrDrift[@]}";do
	echo $nDrift
	sWD="$sOutDir/drift$nDrift"
	mkdir -p $sWD;
	nDrift2=`echo "$nDrift 0.00001" | awk '{printf "%f", $1 + $2}'`
	cat cf.tmpl |  sed 's/%DRIFT1%/'$nDrift'/' | sed 's/%DRIFT2%/'$nDrift2'/' > $sWD/hybrid.cfg
	ln -sf `realpath .`/group1.fa $sWD/
	ln -sf `realpath .`/mixnmatch.sbatch $sWD/
	ln -sf `realpath .`/ancinfer.sbatch $sWD/
	ln -sf `realpath .`/hmm.tmpl $sWD/hmm.cfg
	ln -sf `realpath .`/evaluate.sh $sWD/
done
