#rm Ss*
#cp /fast3/group_crf/home/b20gengjlin3/xipho_simul/*/Sstats.*.txt .
for i in Ss.hm.s.txt;do
	sBase=`basename $i`
	sName=${sBase/%.txt}
	sName=$(echo "$sName" | cut -d'.' -f2)
	echo $sName
	sed -i '/^Gen/!d' $i &&
	sed -i 's/Gen//g' $i &&
#这里要先过滤不为Gen开头的行
	sed -i '/nan/d' $i &&
	awk '$9 >= 0' $i > ${sName}.fixed.txt
	./reg -p ${sName}.fixed.txt -d neutral_data -P 5 -S 4 -b neutral.out -t 0.05 -T 

done
wait
