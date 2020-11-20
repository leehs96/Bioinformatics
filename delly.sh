normal=$1
tumor=$2


delly call -t DEL -q 15 -o ./delly_results/$2.DEL.bcf -g ~/ref/hg37/human_g1k_v37.fasta ./$2.s.md.br.bam ./$1.s.md.br.bam &> $2.DEL.out & 
delly call -t INS -q 15 -o ./delly_results/$2.INS.bcf -g ~/ref/hg37/human_g1k_v37.fasta ./$2.s.md.br.bam ./$1.s.md.br.bam &> $2.INS.out &
delly call -t DUP -q 15 -o ./delly_results/$2.DUP.bcf -g ~/ref/hg37/human_g1k_v37.fasta ./$2.s.md.br.bam ./$1.s.md.br.bam &> $2.DUP.out &
delly call -t INV -q 15 -o ./delly_results/$2.INV.bcf -g ~/ref/hg37/human_g1k_v37.fasta ./$2.s.md.br.bam ./$1.s.md.br.bam &> $2.INV.out &
delly call -t TRA -q 15 -o ./delly_results/$2.TRA.bcf -g ~/ref/hg37/human_g1k_v37.fasta ./$2.s.md.br.bam ./$1.s.md.br.bam &> $2.TRA.out
