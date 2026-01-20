#index
bcftools index *.vcf.gz

#merge
bcftools merge *.vcf.gz -o all.snp.vcf.gz -O z

#index again
bcftools index all.snp.vcf.gz
