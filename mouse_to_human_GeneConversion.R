library(biomaRt)
library(plyr)
setwd('/users/hslee/Pink1')

df <- read.csv('pink1_norm(ensembl).csv')


human = useEnsembl(biomart = "ensembl", dataset="hsapiens_gene_ensembl") 
mouse = useEnsembl(biomart = "ensembl", dataset="mmusculus_gene_ensembl") 


gene <- getLDS(
   attributes = c("ensembl_gene_id"),
   filters = "ensembl_gene_id",
   values = df$X ,
   mart = mouse,
   attributesL = c("mgi_symbol"),
   martL = mouse,
   uniqueRows=T
)


gene <- arrange(gene,Gene.stable.ID)











write.csv(genes,'mouse_to_human.csv',row.names = F)

###################################################################################
#ensembl <- useEnsembl(biomart = "ensembl", dataset="hsapiens_gene_ensembl") 
#chr1genes <- getBM(attributes = c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol',
#                                  'chromosome_name','start_position','end_position'), 
#                   filters ='chromosome_name', values ="1", mart = ensembl)
#
#######################################################################################
#humanx <- unique(genes[, 2])

