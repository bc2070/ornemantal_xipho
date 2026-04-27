dir=$1
codeml=/data/software/paml4.9j/src/codeml
cfg=/data2/projects/jgeng/peak_hm_GO/genespace/Xma/outgroup_prot/cds_outgroups/new_allregion/codeml/codeml.1.ctl
xhel_tre=/data2/projects/jgeng/peak_hm_GO/genespace/Xma/outgroup_prot/cds_outgroups/new_allregion/codeml/tree/${dir}_xhel.tre
xmac_tre=/data2/projects/jgeng/peak_hm_GO/genespace/Xma/outgroup_prot/cds_outgroups/new_allregion/codeml/tree/${dir}_xmac.tre
xvar_tre=/data2/projects/jgeng/peak_hm_GO/genespace/Xma/outgroup_prot/cds_outgroups/new_allregion/codeml/tree/${dir}_xvar.tre

if [ "$dir" == "single_o1" ]; then
    outgroup="gambusia"
elif [ "$dir" == "single_o2" ]; then
    outgroup="poecilia"
else
    outgroup="o1o2" 
fi


for i in /data2/projects/jgeng/peak_hm_GO/genespace/Xma/outgroup_prot/cds_outgroups/new_allregion/${dir}/ENSXMAT*; do
    cd "$i" || continue  
    cp "$cfg" .
    in_1="${i}/${outgroup}_hellerii.codeml.input"
    in_2="${i}/${outgroup}_maculatus.codeml.input"
    in_3="${i}/${outgroup}_variatus.codeml.input"
    out_1="${in_1%.input}.1.out"
    out_2="${in_2%.input}.1.out"
    out_3="${in_3%.input}.1.out"


    sed "s#%INPUT%#$in_1#g" codeml.1.ctl > codeml.xhel.ctl
    sed "s#%TREE%#$xhel_tre#g" codeml.xhel.ctl > codeml.xhel.1.ctl
    sed "s#%OUTPUT%#$out_1#g" codeml.xhel.1.ctl > codeml.xhel.ctl

    sed "s#%INPUT%#$in_2#g" codeml.1.ctl > codeml.xmac.ctl
    sed "s#%TREE%#$xmac_tre#g" codeml.xmac.ctl > codeml.xmac.1.ctl
    sed "s#%OUTPUT%#$out_2#g" codeml.xmac.1.ctl > codeml.xmac.ctl

    sed "s#%INPUT%#$in_3#g" codeml.1.ctl > codeml.xvar.ctl
    sed "s#%TREE%#$xvar_tre#g" codeml.xvar.ctl > codeml.xvar.1.ctl
    sed "s#%OUTPUT%#$out_3#g" codeml.xvar.1.ctl > codeml.xvar.ctl

    rm codeml.x*.1.ctl  
    $codeml codeml.xhel.ctl
    $codeml codeml.xmac.ctl
    $codeml codeml.xvar.ctl
done
wait

