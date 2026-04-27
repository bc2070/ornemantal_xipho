dir=$1


for i in /data2/projects/jgeng/peak_hm_GO/genespace/Xma/outgroup_prot/cds_outgroups/new_allregion/${dir}/*/*.cds;do
	sDir=`dirname $i`
        sBase=`basename $i`
        sName=${sBase/%.cds}
	sOut=${sDir}/${sName}.codeml.input
	echo $i
	echo $sOut
	python3 01.1.codeml.process.py $i > $sOut

done
wait



