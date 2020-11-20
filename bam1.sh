sampleName=$1


echo 'Base recalib'
gatk --java-options -"Xms20G" BaseRecalibrator -R ~/ref/hg37/human_g1k_v37.fasta -I $1.s.md.bam --known-sites ~/ref/hg37/Mills_and_1000G_gold_standard.indels.b37.vcf --known-sites ~/ref/hg37/dbsnp_138.b37.vcf -O $1.s.md.bam.table &&
sleep 5s
gatk --java-options -"Xms20G"  ApplyBQSR  -R ~/ref/hg37/human_g1k_v37.fasta -I $1.s.md.bam --bqsr-recal-file $1.s.md.bam.table -O $1.s.md.br.bam > $sampleName.BR.out 2>&1  
