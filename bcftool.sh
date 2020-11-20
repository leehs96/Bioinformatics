bcf=$1


echo "start SV annotation from bcf combine"                                                                          
bcftools concat -a -O v -o ./results/$1.delly.vcf $1.DEL.bcf $1.DUP.bcf $1.INS.bcf $1.INV.bcf $1.TRA.bcf &> $1.concat.out
echo "done" 
