library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(dplyr)
library(rstatix)
library(RColorBrewer)
library(gplots)
library(ggpubr)


getwd()
setwd("./TCGA")


df <- read.csv('THCA_RNAseq_TCGA.csv', sep = ',', header = T, row.names = 1)
df1 <- read.csv('BRCA_RNAseq_TCGA.csv', sep = ',', header = T, row.names = 1)
df2 <- read.csv('TCGA_SKCM_non_Log.csv', sep = ',', header = T, row.names = 1)
####################################################################################
df_nolog <- read.csv('TCGA_THCA_non_Log.csv', sep = ',', header = T, row.names = 1)

df_ptc_nolog <- select(df_nolog, grep('.01', names(df_nolog)))

df_nt_nolog <- select(df_nolog, grep('.11', names(df_nolog)))



write.csv(df_ptc_nolog, 'PTC_genescoring.csv')
write.csv(df_nt_nolog, 'NT_genescoring.csv')



#####################################################################################

df_row <- df[,1]
rownames(df) <- df_row
df <- df[,-1]


df_m <- melt(as.matrix(df))
names(df_m) <- c("gene","sample","exp")


df_PTC <- select(df, grep('.01', names(df)))

df_cm <- select(df, grep('.06', names(df)))

df_n <- select(df, grep('.11', names(df)))
###############################################################################
df <- read.csv('LUNG_RNAseq_TCGA.csv', sep = ',', header = T, row.names = 1)

df_n <- select(df, grep('.11', names(df)))

df_C <- select(df, grep('.01', names(df)))

df_cp <- df_cp[,-1]
df_ptc_t <- t(df_cp)
df_ptc_t <- as.data.frame(df_ptc_t)

df_nt
df_ptc_t
write.csv(df_ptc_t,'test.csv')

df_high_CD14_nt <- subset(df_ptc_t, CD14 >= 9.87)
df_low_CD14_nt <- subset(df_ptc_t, CD14 < 9.87)

df_high_CD14_nt <- arrange(df_high_CD14_nt, desc(CD14))
df_low_CD14_nt <- arrange(df_low_CD14_nt, CD14)

high_CD14_nt <- t(df_high_CD14_nt)
low_CD14_nt <- t(df_low_CD14_nt)


top10_high_nt <- high_CD14_nt[,c(1:10)]
top10_low_nt <- low_CD14_nt[,c(1:10)]

setwd('../Cibersort/PTC')
write.csv(top10_high_nt, 'top10_high_nt.csv')
write.csv(top10_low_nt, 'top10_low_nt.csv')





###############################################################################

df_ptc <- select(df, grep('.01', names(df)))

df_mptc <- select(df, grep('.06', names(df)))

df_n <- select(df, grep('.11', names(df)))

df_N <-melt(as.matrix(df_n)) ; names(df_N) <- c("gene","sample","exp")
df_PTC <-melt(as.matrix(df_ptc)) ; names(df_PTC) <- c("gene","sample","exp")
df_MPTC <-melt(as.matrix(df_mptc)) ; names(df_MPTC) <- c("gene","sample","exp")

df_N$sample <- 'Normal'
df_PTC$sample <- 'Primary PTC'
df_MPTC$sample <- 'Metastatic PTC'

df_box <- rbind(df_N, df_PTC, df_MPTC)

df_box_CD44 <- filter(df_box, gene == "CD44")

CD44_box <- ggboxplot(
  df_box_CD44,
  x = "sample",
  y = "exp",
  fill = "sample",
  palette = 'jco') + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25))
  ) +
  ylab("Log2(normalized counts + 1)") +
  xlab('') +
  ggtitle("CD44 expression")
CD44_box

stat.test55 <- df_box_CD44 %>% t_test(exp ~ sample)
stat.test55 <- stat.test55 %>% add_xy_position(x = "sample")
CD44_box + stat_pvalue_manual(stat.test55, label = "p.adj.signif", tip.length = 0.01) 
stat.test55

stat.test_ptc <- df_box_CD14 %>%
  t_test(exp ~ sample, ref.group = "Normal") %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test_ptc
stat.test_ptc <- stat.test_ptc %>% add_xy_position(x = "sample")
stat.test_ptc

CD44_box + geom_jitter(alpha = 0.3, color = "#2e4057" , position=position_jitter(0.35)) + stat_pvalue_manual(stat.test55, label = "p.adj.signif", tip.length = 0.01)



df_nonLog <- df

#for(i in 1:nrow(df)) {
  for(j in 1:ncol(df)) {
    
    df_nonLog[j] <- 2^df[j]
  }
#}

write.csv(df_nonLog, 'TCGA_BRCA_non_Log.csv')

ciber <- as.data.frame(t(df))
colnames(ciber) <- ciber[c(1),]
ciber <- ciber[c(-1),]

ciber <- as.data.frame(ciber)

min(ciber$CD14)
max(ciber$CD14)

ciber$CD14 <- as.numeric(ciber$CD14)

df_high_CD14 <- subset(ciber, CD14 >= 9.87)
df_low_CD14 <- subset(ciber, CD14 < 9.87)

df_high_CD14 <- arrange(df_high_CD14, desc(CD14))
df_low_CD14 <- arrange(df_low_CD14, CD14)

high_CD14 <- t(df_high_CD14)
low_CD14 <- t(df_low_CD14)


top10_high <- high_CD14[,c(1:10)]
top10_low <- low_CD14[,c(1:10)]

write.csv(top10_high, 'high_CD14.csv')
write.csv(top10_low, 'low_CD14.csv')



CD14_max <- subset(ciber, CD14 == 13.4237)
CD14_max <- subset(ciber, CD14 == 7.2117)

df_n
df_cp

df_n$X <- rownames(df_n)
df_cp$X <- rownames(df_cp)

df_cp <- df_cp[c(506,1:505)]
df_n <- df_n[c(60,1:59)]

################################################################################################################
getwd()
setwd('./new')
df <- read.csv('THCA_RNAseq_TCGA.csv', sep = ',', header = T)
df <- read.csv('BRCA_RNAseq_TCGA.csv', sep = ',', header = T)
df <- read.csv('SKCM_RNAseq_TCGA.csv', sep = ',', header = T)

df_n <- select(df, grep('.01$', names(df)))
df_n$X <- df$X 
df_n <- df_n[c(506,1:505)]

df_c <- select(df, grep('.06$', names(df)))
df_c$X <- df$X 
df_c <- df_c[c(369,1:368)]

gene_phago <- read.table('phagocytosis_21.txt', header = T)
gene_prog <- read.table('Programmed cell death..txt', header = T)
gene_m1 <- read.table('M1.txt', header = T)
gene_m2 <- read.table('M2.txt', header = T)
gene_tam <- read.table('TAM.txt', header = T)
gene_FCgR_phago <- read.table('FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS.txt', header = T)

M1_df <- merge(df_n, gene_FCgR_phago, by = 1)
M2_df <- merge(df_n, gene_m2, by = 1)
TAM_df <- merge(df_n, gene_tam, by = 1)
apop_df <- merge(df_n, gene_prog, by = 1)
phago_df <- merge(df_n, gene_phago, by = 1)


row.names(apop_df) <- apop_df$X
row.names(phago_df) <- phago_df$X
row.names(TAM_df) <- TAM_df$X
row.names(M1_df) <- M1_df$X
row.names(M2_df) <- M2_df$X


apop_df <- apop_df %>% select(!X)
phago_df <- phago_df %>% select(!X)
TAM_df <- TAM_df %>% select(!X)
M1_df <- M1_df %>% select(!X)
M2_df <- M2_df %>% select(!X)


apop_df_t <- t(apop_df)
phago_df <- t(phago_df)
TAM_df <- t(TAM_df)
M1_df <- t(M1_df)
M2_df <- t(M2_df)


apop_df_t <- as.data.frame(apop_df_t)
phago_df_t <- as.data.frame(phago_df)
TAM_df_t <- as.data.frame(TAM_df)
M1_df_t <- as.data.frame(M1_df)
M2_df_t <- as.data.frame(M2_df)


apop_df_t_d <- arrange(apop_df_t, CD14)
phago_df_t_d <- arrange(phago_df_t, CD14)
TAM_df_t_d <- arrange(TAM_df_t, CD14)
M1_df_t_d <- arrange(M1_df_t, CD14)
M2_df_t_d <- arrange(M2_df_t, CD14)


apop_df_t_d_t <- t(apop_df_t_d)
phago_df_t_d_t <-t(phago_df_t_d)
TAM_df_t_d_t <-t(TAM_df_t_d)
M1_df_t_d_t <-t(M1_df_t_d)
M2_df_t_d_t <-t(M2_df_t_d)

TAM_df_t_d_t_f <- TAM_df_t_d_t[-c(2),]
M1_df_t_d_t_f <- M1_df_t_d_t[-c(1),]
M2_df_t_d_t_f <- M2_df_t_d_t[-c(2),]
phago_df_t_d_t_f <- phago_df_t_d_t[-c(49),]
apop_df_t_d_t_f <- apop_df_t_d_t[-c(4),]

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


library(reshape)
library(pheatmap)
library(dplyr)
library(gplots)

getwd()
setwd('../')

gs <- read.csv('SKCMm.PROJ.csv', sep = ',', header = T)

gs_r <- gs[,1]
rownames(gs) <- gs_r
gs <- gs[,-1]
gst <- t(gs)

gst <- as.data.frame(gst)
gst <- arrange(gst, CD14)

gst <- t(gst)

gst_heat <- gst[-c(6),]



p3 <- pheatmap(
  gst_heat,
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
  legend = F,
)

p3





###






df <- read.csv('THCA_RNAseq_TCGA.csv', sep = ',', header = T)

df_cp <- t(df_cp) 
df_cp <- as.data.frame(df_cp)
df <- t(df)

df_row <- df[,1]
rownames(df) <- df_row
df <- df[,-1]


df <- as.data.frame(df)


df_1 <- df[c('FCGR3A', 'CD36', 'TIMD4', 'MERTK', 'STAB1', 'STAB2', 'ITGB3', 'ITGB5', 'GAS6', 'MFGE8', 'TLR4', 'TLR5', 'TLR6', 'TLR7', 'TLR9', 'AXL', 'RAGE', 'CD300A', 'LRP1'),]
# correlation with CD14

#CD16, CD36, Timd4, MERTK, Stab1, Stab2, Itgb3, Itgb5, Gas6, Mfge8, TLR4, TLR5, TLR6, TLR7, TLR9, AXL, RAGE, CD300a, LPX1, LRP1

#CD16 (FCGR3A)

df_C <- t(df_C)
df_N <- t(df_n)

df_C <- as.data.frame(df_C)
df_N <- as.data.frame(df_N)

CD44 <-ggscatter(df_N, x = "CD44", y = "CD14",
                 add = "reg.line",
                 add.params = list(color = "blue", fill = "gray"),
                 conf.int = T
) 
CD44 + stat_cor(method = 'pearson')



FCGR3A <- ggscatter(df_nt, x = "FCGR3A", y = "CD14",
                  add = "reg.line",
                  add.params = list(color = "blue", fill = "gray"),
                  conf.int = T
) + xlab("CD16")
FCGR3A + stat_cor(method = 'pearson')


setwd('/10.80.51.8/hslee/TCGA/sep/scatter')

OLR1 <- ggscatter(df_cp, x = "OLR1", y = "CD14",
                  add = "reg.line",
                  add.params = list(color = "blue", fill = "gray"),
                  conf.int = T
                  ) 
OLR1 + stat_cor(method = 'pearson')


#for (i in seq_along(list)) {
#  list[i] <- ggscatter(df, x = list[i], y = "CD14",
#                    add = "reg.line",
#                    add.params = list(color = "blue", fill = "gray"),
#                    conf.int = T
#  ) + xlab('CD16')
#  
#  list[i] <- list[i] + stat_cor(method = 'pearson', label.x = 8, label.y = 13)
#}

#CD36

CD36 <- ggscatter(df_cp, x = "CD36", y = "CD14",
                  add = "reg.line",
                  add.params = list(color = "blue", fill = "gray"),
                  conf.int = T
) #+ xlab('CD36')
CD36 + stat_cor(method = 'pearson')


#Timd4

TIMD4 <- ggscatter(df_cp, x = "TIMD4", y = "CD14",
                  add = "reg.line",
                  add.params = list(color = "blue", fill = "gray"),
                  conf.int = T
) + xlab('TIMD4')
TIMD4 + stat_cor(method = 'pearson')


#MERTK

MERTK <- ggscatter(df_cp, x = "MERTK", y = "CD14",
                  add = "reg.line",
                  add.params = list(color = "blue", fill = "gray"),
                  conf.int = T
) + xlab('MERTK')
MERTK + stat_cor(method = 'pearson')


#Stab1

STAB1 <- ggscatter(df_cp, x = "STAB1", y = "CD14",
                  add = "reg.line",
                  add.params = list(color = "blue", fill = "gray"),
                  conf.int = T
) + xlab('STAB1')
STAB1 + stat_cor(method = 'pearson')


#Stab2

STAB2 <- ggscatter(df_cp, x = "STAB2", y = "CD14",
                  add = "reg.line",
                  add.params = list(color = "blue", fill = "gray"),
                  conf.int = T
) + xlab('STAB2')
STAB2 + stat_cor(method = 'pearson')


#Itgb3

ITGB3 <- ggscatter(df_cp, x = "ITGB3", y = "CD14",
                  add = "reg.line",
                  add.params = list(color = "blue", fill = "gray"),
                  conf.int = T
) + xlab('ITGB3')
ITGB3 + stat_cor(method = 'pearson')

#Itgb5

ITGB5 <- ggscatter(df_cp, x = "ITGB5", y = "CD14",
                  add = "reg.line",
                  add.params = list(color = "blue", fill = "gray"),
                  conf.int = T
) + xlab('ITGB5')
ITGB5 + stat_cor(method = 'pearson')

#Gas6

GAS6 <- ggscatter(df_cp, x = "GAS6", y = "CD14",
                  add = "reg.line",
                  add.params = list(color = "blue", fill = "gray"),
                  conf.int = T
) + xlab('GAS6')
GAS6 + stat_cor(method = 'pearson')

#Mfge8

MFGE8 <- ggscatter(df_cp, x = "MFGE8", y = "CD14",
                  add = "reg.line",
                  add.params = list(color = "blue", fill = "gray"),
                  conf.int = T
) + xlab('MFGE8')
MFGE8 + stat_cor(method = 'pearson')

#TLR4

TLR4 <- ggscatter(df_cp, x = "TLR4", y = "CD14",
                  add = "reg.line",
                  add.params = list(color = "blue", fill = "gray"),
                  conf.int = T
) + xlab('TLR4')
TLR4 + stat_cor(method = 'pearson')

#TLR5

TLR5 <- ggscatter(df_cp, x = "TLR5", y = "CD14",
                  add = "reg.line",
                  add.params = list(color = "blue", fill = "gray"),
                  conf.int = T
)
TLR5 + stat_cor(method = 'pearson')


#TLR6

TLR6 <- ggscatter(df_cp, x = "TLR6", y = "CD14",
                   add = "reg.line",
                   add.params = list(color = "blue", fill = "gray"),
                   conf.int = T
)
TLR6 + stat_cor(method = 'pearson')

#TLR7

TLR7 <- ggscatter(df_cp, x = "TLR7", y = "CD14",
                  add = "reg.line",
                  add.params = list(color = "blue", fill = "gray"),
                  conf.int = T
)
TLR7 + stat_cor(method = 'pearson')
#TLR9

TLR9 <- ggscatter(df_cp, x = "TLR9", y = "CD14",
                  add = "reg.line",
                  add.params = list(color = "blue", fill = "gray"),
                  conf.int = T
)
TLR9 + stat_cor(method = 'pearson')

#AXL

AXL <- ggscatter(df_cp, x = "AXL", y = "CD14",
                  add = "reg.line",
                  add.params = list(color = "blue", fill = "gray"),
                  conf.int = T
)
AXL + stat_cor(method = 'pearson')

#RAGE

RAGE <- ggscatter(df_cp, x = "RAGE", y = "CD14",
                  add = "reg.line",
                  add.params = list(color = "blue", fill = "gray"),
                  conf.int = T
)
RAGE + stat_cor(method = 'pearson')
#CD300A

CD300A <- ggscatter(df_cp, x = "CD300A", y = "CD14",
                  add = "reg.line",
                  add.params = list(color = "blue", fill = "gray"),
                  conf.int = T
)
CD300A + stat_cor(method = 'pearson')
#LPX1

CD47 <- ggscatter(df_cp, x = "CD47", y = "CD14",
                  add = "reg.line",
                  add.params = list(color = "blue", fill = "gray"),
                  conf.int = T
)
CD47 + stat_cor(method = 'pearson')

#LRP1

LRP1 <- ggscatter(df_cp, x = "LRP1", y = "CD14",
                  add = "reg.line",
                  add.params = list(color = "blue", fill = "gray"),
                  conf.int = T
)
LRP1 + stat_cor(method = 'pearson')


df <- read.csv('TCGA_THCA_non_Log.csv', sep = ',', header = T, row.names = 1)
df1 <- read.csv('TCGA_BRCA_non_Log.csv', sep = ',', header = T, row.names = 1)
df2 <- read.csv('TCGA_SKCM_non_Log.csv', sep = ',', header = T, row.names = 1)


PTC_t <- select(df, grep('.01', names(df)))
BRCA_t <- select(df1, grep('.01', names(df1)))
SKCM_t <- select(df2, grep('.01', names(df2)))
SKCM_m <- select(df2, grep('.06', names(df2)))

write.csv(PTC_t,'PTC_t.csv')
write.csv(BRCA_t,'BRCA_t.csv')
write.csv(SKCM_t,'SKCM_t.csv')
write.csv(SKCM_m,'SKCM_m.csv')
