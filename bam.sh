sampleName=$1

echo 'bwa'
bwa mem -M -t 6 -R "@RG\tID:$1\tLB:1\tSM:$1\tPL:ILLUMINA" ~/ref/hg37/human_g1k_v37.fasta.gz ~/../data/CancerEvolution_WGS/fastq/$1_DNA_R1.fastq.gz ~/../data/CancerEvolution_WGS/fastq/$1_DNA_R2.fastq.gz > $1.sam 2> $1.bwamem.out &&
echo 'sam to bam'
samtools view -Sb -@ 6 -o $1.bam $1.sam > $1.StoB.out 2>&1 &&
rm $1.sam
echo 'sort'
samtools sort -@ 6 -o $1.s.bam $1.bam > $1.sort.out 2>&1 &&
rm $1.bam
echo 'Markdup' &&
java -XX:ParallelGCThreads=6 -XX:ConcGCThreads=6  -Xms12g -Xmx20g -jar ~/anaconda3/pkgs/picard-2.23.8-0/share/picard-2.23.8-0/picard.jar MarkDuplicates -REMOVE_DUPLICATES true -REMOVE_SEQUENCING_DUPLICATES true -I $1.s.bam -O $1.s.md.bam -M $1.s.metrics.txt -VALIDATION_STRINGENCY LENIENT -CREATE_INDEX true &&
echo 'Base recalib'
gatk --java-options -"Xms20G" BaseRecalibrator -R ~/ref/hg37/human_g1k_v37.fasta.gz -I $1.s.md.bam --known-sites ~/ref/hg37/Mills_and_1000G_gold_standard.indels.b37.vcf --known-sites ~/ref/hg37/dbsnp_138.b37.vcf -O $1.s.md.bam.table 
gatk --java-options -"Xms20G"  ApplyBQSR  -R ~/ref/hg37/human_g1k_v37.fasta.gz -I $1.s.md.bam --bqsr-recal-file $1.s.md.bam.table -O $1.s.md.br.bam > $sampleName.BR.out 2>&1
rm $1.s.md.bam
rm $1.s.md.bam.table
