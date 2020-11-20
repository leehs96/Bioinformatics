normal=$1
tumor=$2


java -Xms20G -jar ~/bin/varscan-2.4.4.jar somatic ./mpileup/$1.mpileup ./mpileup/$2.mpileup ./varscan_results/$2.varscan --min-var-freq 0.01 --output-vcf 1 &> ./varscan_results/$2_varscan.out &&
mv ./varscan_results/$2_varscan.out ./varscan_results/$2.varscan.success
