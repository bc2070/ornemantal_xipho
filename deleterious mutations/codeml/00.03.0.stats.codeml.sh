dir=$1


for i in /data2/projects/jgeng/peak_hm_GO/genespace/Xma/outgroup_prot/cds_outgroups/new_allregion/${dir}/*/*codeml.1.out;do
	sDir=`dirname $i`
	sBase=`basename $i`
        sName=${sBase/%.codeml.1.out}
	sSp=${sName#*_}
	sID=${sDir##*/}
#	echo $sSp
#	echo $sID

	value=$(awk 'NR==89 {print $11}' $i)	#dN
#	value=$(awk 'NR==89 {print $8}' $i)	#dN/dS omega
#pairwise cal
#	value=$(awk '/w ratios as labels for TreeView:/ {getline; if (match($0, /Xiphophorus_[a-zA-Z_]+ #*([0-9.]+)/, arr)) print arr[1]; else print ""}' "$i")
#user tree cal

#	echo $value
	echo -e "$sID\t$sSp\t$value" >> ${dir}.stats.txt

done
wait

