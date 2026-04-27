gff=/data2/projects/jgeng/peak_hm_GO/genespace/Xma/Xiphophorus_maculatus.X_maculatus-5.0-male.110.gff3
fasta=/data2/projects/jgeng/peak_hm_GO/genespace/Xma/Xiphophorus_variatus_pseudogenome2.fasta
name=Xiphophorus_variatus
/data/software/PASApipeline.v2.4.1/misc_utilities/gff3_file_single_longest_isoform.pl $gff > ${name}.longest.isoform.gff3
/data/software/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl ${name}.longest.isoform.gff3 $fasta cDNA > ${name}.longest.isoform.mrna.gff3
/data/software/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl ${name}.longest.isoform.gff3 $fasta prot > ${name}.longest.isoform.prot.gff3
/data/software/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl ${name}.longest.isoform.gff3 $fasta CDS > ${name}.longest.isoform.cds.gff3
