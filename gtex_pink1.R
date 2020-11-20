setwd('/users/hslee/gtex')

library(dplyr)
library(ggplot2)
library(reshape)
library(ggpubr)
library(rstatix)
library(dplyr)

#raw data

gtex_tpm <- read.table('GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.txt',sep = '\t',header = T)
gtex_anno <- read.csv('gtex_anno.csv', sep = ",", header = T)
gtex_pink1 <- read.csv('gtex_tpm_pink1.csv', header = T, row.names = 1)
#parsing
#gtex_tpm_pink1 <- gtex_tpm %>% filter(gtex_tpm$Description == "PINK1")
gtex_anno_thyroid <- gtex_anno %>% filter(gtex_anno$SMTS == 'Thyroid')

test <-gtex_anno %>% filter(gtex_anno$SMTS == 'Thyroid')

gtex_anno_thyroid$SMPTHNTS <- ifelse(grepl('thyroiditis', gtex_anno_thyroid$SMPTHNTS), 'LT','No-LT')




#write.csv(gtex_tpm_pink1,file = "gtex_tpm_pink1.csv", row.names = 1)
melt_gtex_pink1 <- melt(as.matrix(gtex_pink1))
melt_gtex_pink1 <- na.omit(melt_gtex_pink1)

names(melt_gtex_pink1) <- c('PINK1','SAMPID','TPM')
gtex_anno_thyroid$SAMPID <- gsub('-','.',gtex_anno_thyroid$SAMPID)

gtex_data <- merge(melt_gtex_pink1, gtex_anno_thyroid, by = 'SAMPID')
gtex_data <- gtex_data %>% select(SAMPID, PINK1, SMPTHNTS,TPM)

#----------------------------------------------------------------------

#box plot 

gtex_box <- ggboxplot(
  gtex_data,
  x = "SMPTHNTS",
  y = "TPM",
  fill = "SMPTHNTS",
  palette = c("#00798c","#d1495b"),
  notch = TRUE) + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25))
  ) +
  ylab("normalized counts") +
  xlab("GTEx Analysis  T-test, p = 0.0000000419") +
  ggtitle("PINK1 expression")
gtex_box

#stat test

gtex_stat <- gtex_data %>%
  t_test(TPM ~ SMPTHNTS) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
gtex_stat


gtex_stat <- gtex_stat %>% add_xy_position(x = 'SMPTHNTS')


#star
gtex_box +
  stat_pvalue_manual(
    gtex_stat,
    label = "p.adj.signif",
    vjust = -0.5,
    bracket.nudge.y = 2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#t-test
gtex_box +
  stat_pvalue_manual(
    gtex_stat,
    label = "p.adj.signif",
    vjust = -0.5,
    bracket.nudge.y = 2
  ) +
  geom_jitter(alpha = 0.4, color = "#2e4057" , position=position_jitter(0.35)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) 
