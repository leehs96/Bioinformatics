# library

library(pheatmap)
library(TCGAbiolinks)
library(heatmaply)
library(RColorBrewer)
library(gplots)
library(ggpubr)
library(viridis)

## find dataset
# 1. Search data in GDC
query <- GDCquery(project="TCGA-THCA",
                  data.category="Transcriptome Profiling",
                  data.type="Gene Expression Quantification",
                  workflow.type="HTSeq - FPKM-UQ",
                  )

# 2. Download from GDC repository

GDCdownload(query)

# 3. Make R object from the downloaded data

data <- GDCprepare(query)

# 4. Extract Gene expression matrix

library(SummarizedExperiment)
eset <- assay(data)

# 5. Save the matrix as .csv format

write.csv(eset, file="THCA_GE.csv")

# 안쓸래 이거 <- Gene symbol 데이터 없음 다 ensembl ID 로 되어있어서 변환해줘야함
###############################################################################################################################################

getwd()
setwd('~/TAM')

df <- read.csv('BRCA_RNAseq_TCGA.csv', sep = ',', header = T)
#df <- read.csv('BRCA_RNAseq_TCGA.csv', sep = ',', header = T)
#df <- read.csv('SKCM_RNAseq_TCGA.csv', sep = ',', header = T)

df_n <- select(df, grep('.01$', names(df)))
df_n$X <- df$X 
df_n <- df_n[c(1098,1:1097)]

df_c <- select(df, grep('.06$', names(df)))
df_c$X <- df$X 
df_c <- df_c[c(369,1:368)]

M2 <- read.table('M2_gene_list.txt', header = T)
gene_MMP <- read.table('Programmed cell death..txt', header = T)
gene_ADAM <- read.table('M1.txt', header = T)
gene_Succinate <- read.table('M2.txt', header = T)
gene_Ala_Asp <- read.table('TAM.txt', header = T)
gene_chemoK. <- read.table('FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS.txt', header = T)
gene_effero <- read.table('FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS.txt', header = T)

standard = 'CD163'

Heat <- function(df, gene_set, standard){
  gene_set[nrow(gene_set) + 1,] = standard
  gene_df <- merge(df, gene_set, by = 1)
  row.names(gene_df) <- gene_df$X
  gene_df <- gene_df %>% select(!X)
  gene_df_t <- t(gene_df)
  gene_df_t <- as.data.frame(gene_df_t)
  gene_df_t_d <- gene_df_t %>% arrange(gene_df_t[names(gene_df_t)==standard])
  gene_df_t_d_t <- t(gene_df_t_d)
  gene_df_t_d_t_f <- gene_df_t_d_t[!rownames(gene_df_t_d_t)==standard,]
  
  p <- pheatmap(
    gene_df_t_d_t_f,
    color = greenred(11),
    breaks = c(-1,-0.9,-0.7,-0.5,-0.2,0,0.2,0.5,0.7,0.9,1),
    clustering_method = "ward.D2",
    border_color = NA,
    clustering_distance_rows = "correlation",
    fontsize_row = 6,
    cluster_rows = T,
    cluster_cols = F,
    scale = "row",
    show_colnames = F,
    show_rownames = T,
    legend = T
  )
  
}

Heat(df_n,M2,standard)




M2_df <- merge(df_n, M2, by = 1)


row.names(M2_df) <- M2_df$X



M2_df <- M2_df %>% select(!X)



M2_df_t <- t(M2_df)


M2_df_t <- as.data.frame(M2_df_t)



M2_df_t_d <- arrange(M2_df_t, M2_df_t[names(M2_df_t)=='ARSK'])

head(names(M2_df_t))

apop_df_t_d_t <- t(apop_df_t_d)


M2_df_t_d_t_f <- TAM_df_t_d_t[-rownames(TAM_df_t_d_t)=='CD163',]
M1_df_t_d_t_f <- M1_df_t_d_t[-c(1),]


p <- pheatmap(
  M1_df_t_d_t,
  color = greenred(11),
  breaks = c(-1,-0.9,-0.7,-0.5,-0.2,0,0.2,0.5,0.7,0.9,1),
  clustering_method = "ward.D2",
  border_color = NA,
  clustering_distance_rows = "correlation",
  fontsize_row = 6,
  cluster_rows = T,
  cluster_cols = F,
  scale = "row",
  show_colnames = F,
  show_rownames = T,
  legend = T
)
p

pr <- rownames(phago_df_t_d_t_f[p$tree_row[["order"]],])
pr_d<-as.data.frame(pr)
write.table(pr_d, 'SKCM_metastatic_phago.txt', quote = F, row.names = F)



