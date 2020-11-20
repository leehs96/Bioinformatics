library(RColorBrewer)
library(gplots)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(ggpubr)
library(reshape)

## Gene set
getwd()
setwd('~/Gene_set')

#########################################################################
path <- "/users/hslee/Gene_set/"
data1 <-list.files(path = path, pattern = '*.txt', full.names = T)
############################################
for (i in seq_along(data1)){
  a <- gsub(path,"",data1[i])
  b <- gsub('.txt','',a)
  name <- gsub('/','',b)
  assign(name, read.table(data1[i], header = T))
}
############################################################
geneset <- read.csv('Immune_2020.11.17.csv', sep = ',', header = T)

############################################
for (i in seq_along(geneset)){
  name <- (colnames(geneset)[i])
  assign(name, as.data.frame(geneset[[i]][geneset[[i]] != '']))
}
############################################
path <- "/users/hslee/TCGA_RNA_data/oh/"
data2 <-list.files(path = path, pattern = '*.csv', full.names = T)
############################################
for (i in seq_along(data2)){
  a <- gsub(path,"",data2[i])
  b <- gsub('/','',a)
  name <- (gsub("_RNAseq_TCGA.csv","_M",b))
  a <- read.csv(data2[i], header = T , sep = ',')
  gene <- a[,1] 
  a <- dplyr::select(a, grep('.01$', names(a)))
  a$X <- gene
  a <- a %>% relocate(X)
  as.data.frame(assign(name, a))
}
############################################
Heat <- function(df, gene_set, standard){
  ifelse(!(standard %in% gene_set[[1]]),gene_set[nrow(gene_set) + 1,] <- standard,NA)
  gene_df <- merge(df, gene_set, by= 1 )
  row.names(gene_df) <- gene_df$X
  gene_df <- gene_df %>% select(!X)
  gene_df_t <- t(gene_df)
  gene_df_t <- as.data.frame(gene_df_t)
  gene_df_t_d <- gene_df_t %>% arrange(gene_df_t[names(gene_df_t)==standard])
  gene_df_t_d_t <- t(gene_df_t_d)
  gene_df_t_d_t_f <- gene_df_t_d_t[!rownames(gene_df_t_d_t)==standard,]
  
  p <- pheatmap(
    gene_df_t_d_t_f,
    color = greenred(5),
    breaks = c(-1,-0.6,0,0.6,1),
    clustering_method = "ward.D2",
    border_color = NA,
    clustering_distance_rows = "correlation",
    fontsize_row = 6,
    cluster_rows = T,
    cluster_cols = F,
    scale = "row",
    show_colnames = F,
    show_rownames = T,
    legend = F,
  )
}
#list#######################################
cancer <- list(BRCA,THCA,SKCM,SKCM_M,LUNG,KIRP,KIRC)
names(cancer) <-c('BRCA','THCA','SKCM','SKCM_M', 'LUNG', 'KIRP' , 'KIRC')
flag <- c("CD14")
geneset <- list(Th1,B_cells,NK_cells,CD8_T_cells)
names(geneset) <-c('Th1','B_cells','NK_cells','CD8_T_cells')
############################################
setwd('~/오병철_교수님/immune_heat/CD14/acute_inf')

for (i in 1:length(cancer)){
  for (j in 1:length(geneset)){
    png(filename=paste0(names(cancer)[[i]],'_',names(geneset)[[j]],'_',flag,'.png'),width=850,height=750,units = "px", bg="white" , res = 120,)
    Heat(cancer[[i]],geneset[[j]],flag)
    dev.off() 
  }
}




