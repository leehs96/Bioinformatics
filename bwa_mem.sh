echo "start"
bwa mem -M -t 5 -R "@RG\tID:ADB_Blood\tLB:1\tSM:ADB_Blood\tPL:ILLUMINA" ./ref/hg38_v0_Homo_sapiens_assembly38.fasta ../data/CancerEvolution_WGS/fastq/ADB_Blood_DNA_R1.fastq.gz ../data/CancerEvolution_WGS/fastq/ADB_Blood_DNA_R2.fastq.gz > ADB_Blood.sam 2> ADB_Blood.bwamem.out 
bwa mem -M -t 5 -R "@RG\tID:ADB_LL\tLB:1\tSM:ADB_LL\tPL:ILLUMINA" ./ref/hg38_v0_Homo_sapiens_assembly38.fasta ../data/CancerEvolution_WGS/fastq/ADB_LL_DNA_R1.fastq.gz ../data/CancerEvolution_WGS/fastq/ADB_LL_DNA_R2.fastq.gz > ADB_LL.sam 2> ADB_LL.bwamem.out
echo "success...........1/14"
bwa mem -M -t 5 -R "@RG\tID:ADB_LM\tLB:1\tSM:ADB_LM\tPL:ILLUMINA" ./ref/hg38_v0_Homo_sapiens_assembly38.fasta ../data/CancerEvolution_WGS/fastq/ADB_LM_DNA_R1.fastq.gz ../data/CancerEvolution_WGS/fastq/ADB_LM_DNA_R2.fastq.gz  > ADB_LM.sam 2> ADB_LM.bwamem.out
echo "success...........2/14"
bwa mem -M -t 5 -R "@RG\tID:ADB_UL\tLB:1\tSM:ADB_UL\tPL:ILLUMINA" ./ref/hg38_v0_Homo_sapiens_assembly38.fasta ../data/CancerEvolution_WGS/fastq/ADB_UL_DNA_R1.fastq.gz ../data/CancerEvolution_WGS/fastq/ADB_UL_DNA_R2.fastq.gz > ADB_UL.sam 2> ADB_UL.bwamem.out
echo "success...........3/14"
bwa mem -M -t 5 -R "@RG\tID:ADB_UM\tLB:1\tSM:ADB_UM\tPL:ILLUMINA" ./ref/hg38_v0_Homo_sapiens_assembly38.fasta ../data/CancerEvolution_WGS/fastq/ADB_UM_DNA_R1.fastq.gz ../data/CancerEvolution_WGS/fastq/ADB_UM_DNA_R2.fastq.gz > ADB_UM.sam 2> ADB_UM.bwamem.out
echo "success...........4/14"
bwa mem -M -t 5 -R "@RG\tID:LMS_Blood\tLB:1\tSM:LMS_Blood\tPL:ILLUMINA" ./ref/hg38_v0_Homo_sapiens_assembly38.fasta ../data/CancerEvolution_WGS/fastq/LMS_Blood_DNA_R1.fastq.gz ../data/CancerEvolution_WGS/fastq/LMS_Blood_DNA_R2.fastq.gz > LMS_Blood.sam 2> LMS_Blood.bwamem.out
echo "success...........5/14"
bwa mem -M -t 5 -R "@RG\tID:LMS_Cancer\tLB:1\tSM:LMS_Cancer\tPL:ILLUMINA" ./ref/hg38_v0_Homo_sapiens_assembly38.fasta ../data/CancerEvolution_WGS/fastq/LMS_Cancer_DNA_R1.fastq.gz ../data/CancerEvolution_WGS/fastq/LMS_Cancer_DNA_R2.fastq.gz > LMS_Cancer.sam 2> LMS_Cancer.bwamem.out 
echo "success...........6/14"
bwa mem -M -t 5 -R "@RG\tID:LMS_LN_3R1\tLB:1\tSM:LMS_LN_3R1\tPL:ILLUMINA" ./ref/hg38_v0_Homo_sapiens_assembly38.fasta ../data/CancerEvolution_WGS/fastq/LMS_LN_3R1_DNA_R1.fastq.gz ../data/CancerEvolution_WGS/fastq/LMS_LN_3R1_DNA_R2.fastq.gz > LMS_LN_3R1.sam 2> LMS_LN_3R1.bwamem.out 
echo "success...........7/14"
bwa mem -M -t 5 -R "@RG\tID:LMS_LN_3R1\tLB:1\tSM:LMS_LN_3L1\tPL:ILLUMINA" ./ref/hg38_v0_Homo_sapiens_assembly38.fasta ../data/CancerEvolution_WGS/fastq/LMS_LN_3L1_DNA_R1.fastq.gz ../data/CancerEvolution_WGS/fastq/LMS_LN_3L1_DNA_R2.fastq.gz > LMS_LN_3L1.sam 2> LMS_LN_3L1.bwamem.out 
echo "success...........8/14"
bwa mem -M -t 5 -R "@RG\tID:LMS_Recur\tLB:1\tSM:LMS_Recur\tPL:ILLUMINA" ./ref/hg38_v0_Homo_sapiens_assembly38.fasta ../data/CancerEvolution_WGS/fastq/LMS_Recur_DNA_R1.fastq.gz ../data/CancerEvolution_WGS/fastq/LMS_Recur_DNA_R2.fastq.gz > LMS_Recurr.sam 2> LMS_Recurr.bwamem.out 
echo "success...........9/14"
bwa mem -M -t 5 -R "@RG\tID:SJH_Blood\tLB:1\tSM:SJH_Blood\tPL:ILLUMINA" ./ref/hg38_v0_Homo_sapiens_assembly38.fasta ../data/CancerEvolution_WGS/fastq/SJH_Blood_DNA_R1.fastq.gz ../data/CancerEvolution_WGS/fastq/SJH_Blood_DNA_R2.fastq.gz > SJH_Blood.sam 2> SJH_Blood.bwamem.out 
echo "success...........10/14"
bwa mem -M -t 5 -R "@RG\tID:SJH_Cen1\tLB:1\tSM:SJH_Cen1\tPL:ILLUMINA" ./ref/hg38_v0_Homo_sapiens_assembly38.fasta ../data/CancerEvolution_WGS/fastq/SJH_Cen1_DNA_R1.fastq.gz ../data/CancerEvolution_WGS/fastq/SJH_Cen1_DNA_R2.fastq.gz > SJH_Cen1.sam 2> SJH_Cen1.bwamem.out  
echo "success...........11/14"
bwa mem -M -t 5 -R "@RG\tID:SJH_Isth\tLB:1\tSM:SJH_Isth\tPL:ILLUMINA" ./ref/hg38_v0_Homo_sapiens_assembly38.fasta ../data/CancerEvolution_WGS/fastq/SJH_Isth_DNA_R1.fastq.gz ../data/CancerEvolution_WGS/fastq/SJH_Isth_DNA_R2.fastq.gz > SJH_Isth.sam 2> SJH_Isth.bwamem.out  
echo "success...........12/14"
bwa mem -M -t 5 -R "@RG\tID:SJH_Lat\tLB:1\tSM:SJH_Lat\tPL:ILLUMINA" ./ref/hg38_v0_Homo_sapiens_assembly38.fasta ../data/CancerEvolution_WGS/fastq/SJH_Lat_DNA_R1.fastq.gz ../data/CancerEvolution_WGS/fastq/SJH_Lat_DNA_R2.fastq.gz > SJH_Lat.sam 2> SJH_Lat.bwamem.out  
echo "success...........13/14"
bwa mem -M -t 5 -R "@RG\tID:SJH_LN\tLB:1\tSM:SJH_LN\tPL:ILLUMINA" ./ref/hg38_v0_Homo_sapiens_assembly38.fasta ../data/CancerEvolution_WGS/fastq/SJH_LN_DNA_R1.fastq.gz ../data/CancerEvolution_WGS/fastq/SJH_LN_DNA_R2.fastq.gz > SJH_LN.sam 2> SJH_LN.bwamem.out  
echo "success...........14/14"
echo "all done"




