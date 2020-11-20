#!/usr/bin/env bash

total_files=`find ../data/CancerEvolutionWGS/fastq -name '*.fastq.gz' | wc -l`
arr=( $(ls ../data/CancerEvolutionWGS/fastq/*_*_*.fastq.gz) )
echo "mapping started" >> map.log
echo "---------------" >> map.log

for((i=0; i<$total_files; i+=2))
	{
		sample_name=`echo ${arr[$i]} | awk -F "_" '{print $1}'`
		echo "[mapping running for] $sample_name"
		printf "\n"
		echo "bwa mem -t 6 ref/hg38_v0_Homo_sapiens_assembly38.fasta ${arr[$i]} ${arr[$i+1]} > $sample_name.sam" >> map.log
		bwa mem -t 6 ref/hg38_v0_Homo_sapiens_assembly38.fasta ${arr[$i]} ${arr[$i+1]} > $sample_name.sam
	}
