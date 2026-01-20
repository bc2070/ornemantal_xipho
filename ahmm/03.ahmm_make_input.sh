py=02.ahmm_vcf2ahmm_xiphophorus.py

map=APS2population_mappings

gzip -d all.snp.vcf.gz

wait

python3 $py -v all.snp.vcf  -s $map -g 0 -r 6.67e-8 > APS_all_aim_input

