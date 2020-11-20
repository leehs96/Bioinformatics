normal=$1
tumor=$2

set -e

~/anaconda3/pkgs/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam ../$1.s.md.br.bam --tumorBam ../$2.s.md.br.bam --ref ~/ref/hg37/human_g1k_v37.fasta --runDir $2 &> $2.strelka2.out

sleep 5s


/users/hslee/cancer/WGS_hg37/strelka_results/$2/runWorkflow.py  -m local -j 9 &> $2.strelka2.runWorkflow.out

