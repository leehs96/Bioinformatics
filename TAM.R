library(RColorBrewer)
library(gplots)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(ggpubr)
library(reshape)
library(data.table)
library(utils)
library(tidy)
library(circlize)

## Gene set
getwd()
setwd('~/TAM/geneset')


succinate_GOI <- read.table('succci_goi.txt', header = T) 
ADAM <- read.table('ADAM.txt',header = T)
MMP <- read.table('MMP.txt',header = T)
chemo <- read.table('chemokine.txt',header = T)
effero <- read.table('efferocytosis.txt',header = T)
G_S_T <- read.table('glycine_serine_threonine_metabolism.txt',header = T)
A_D_E <- read.table('KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM.txt',header = T)
M2 <- read.table('M2_gene_list.txt',header = T)

## Dataset & criterion
setwd('~/TCGA_RNA_data')

BRCA1 <- read.csv('BRCA_RNAseq_TCGA.csv',header = T)
LUNG1 <- read.csv('LUNG_RNAseq_TCGA.csv',header = T)
PRAD1 <- read.csv('PRAD_RNAseq_TCGA.csv',header = T)
THCA1 <- read.csv('THCA_RNAseq_TCGA.csv',header = T)
UCEC1 <- read.csv('UCEC_RNAseq_TCGA.csv',header = T)
OV1 <- read.csv('OV_RNAseq_TCGA.csv',header = T)

Gene <- BRCA1$X

## only primary cancer

BRCA <- dplyr::select(BRCA1, grep('.01$', names(BRCA1)))
LUNG <- select(LUNG1, grep('.01$', names(LUNG1)))
PRAD <- select(PRAD1, grep('.01$', names(PRAD1)))
THCA <- select(THCA1, grep('.01$', names(THCA1)))
UCEC <- select(UCEC1, grep('.01$', names(UCEC1)))
OV <- select(OV1, grep('.01$', names(OV1)))

BRCA$X <- Gene
LUNG$X <- Gene
PRAD$X <- Gene
THCA$X <- Gene
UCEC$X <- Gene
OV$X <- Gene

BRCA <- BRCA %>% relocate(X)
LUNG <- LUNG %>% relocate(X)
PRAD <- PRAD %>% relocate(X)
THCA <- THCA %>% relocate(X)
UCEC <- UCEC %>% relocate(X)
OV <- OV %>% relocate(X)






## Heatmap function - one reference
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
    color = colorRampPalette(c('deepskyblue3','grey20','firebrick1'))(11),
    #color = greenred(11),
    breaks = c(-1,-0.8,-0.6,-0.3,-0.2,0,0.2,0.3,0.6,0.8,1),
    clustering_method = "ward.D2",
    border_color = NA,
    clustering_distance_rows = "correlation",
    fontsize_row = 6,
    cluster_rows = T,
    cluster_cols = F,
    scale = "row",
    show_colnames = F,
    show_rownames = T,
    legend = F
  )
  
}

cancer <- list(BRCA,THCA,LUNG,OV,PRAD,UCEC)
names(cancer) <-c('BRCA','THCA','LUNG','OV','PRAD','UCEC')
flag <- c('CD68')
geneset <- list(ADAM, MMP, G_S_T, M2, A_D_E, effero, chemo, succinate_GOI)
names(geneset) <-c('ADAM', 'MMP', 'G_S_T', 'M2', 'A_D_E', 'effero', 'chemo', 'succinate_GOI')

setwd('~/TAM/CD68')

for (i in 1:length(cancer)){
  for (j in 1:length(geneset)){
    png(filename=paste0(names(cancer)[[i]],'_',names(geneset)[[j]],'_',flag,'.png'),width=850,height=750,units = "px", bg="white" , res = 120,)
    Heat(cancer[[i]],geneset[[j]],flag)
    dev.off() 
    
  }
}
######################################################
#gene scoring
setwd('~/TAM/genescoring/out')

path <- '/users/hslee/TAM/genescoring/out'
gs <-list.files(path = path, pattern = '*.csv', full.names = T)

#####################################################################
for (i in seq_along(gs)){
  a <- gsub(path,"",gs[i])
  b <- gsub('/','',a)
  name <- (gsub(".csv","_ssgsea",b))
  assign(name, read.csv( gs[i], header = T))
}#######################################

cancer <- list(BRCA,THCA,LUNG,OV,PRAD,UCEC)
names(cancer) <-c('BRCA','THCA','LUNG','OV','PRAD','UCEC')



for (i in 1:length(cancer)){
  paste0(names(cancer)[[1]],'_gs')
  <- filter(cancer[[i]], X == 'CD163')
}


BRCA_gs <- filter(BRCA, X == 'CXCL16')
THCA_gs <- filter(THCA, X == 'CXCL16')
LUNG_gs <- filter(LUNG, X == 'CXCL16')
OV_gs <- filter(OV, X == 'CXCL16')
PRAD_gs <- filter(PRAD, X == 'CXCL16')
UCEC_gs <- filter(UCEC, X == 'CXCL16')


BRCA_gs <- BRCA_gs %>% relocate(X)
THCA_gs <- THCA_gs %>% relocate(X)
LUNG_gs <- LUNG_gs %>% relocate(X)
OV_gs <- OV_gs %>% relocate(X)
PRAD_gs <- PRAD_gs %>% relocate(X)
UCEC_gs <- UCEC_gs %>% relocate(X)




Heat_gs <- function(df_gs_ssgsea, df_gs, standard){
 
  gs <- rbind(df_gs_ssgsea, df_gs)
  row <- gs$X
  rownames(gs) <- row
  gs <- gs[,-1]
  
  gst <- t(gs)
  
  gst <- as.data.frame(gst)
  gst <- gst %>% arrange(gst[names(gst)==standard])
  
  gst <- t(gst)
  
  gst_heat <- gst[!rownames(gst)==standard,]
  
  
  
  p3 <- pheatmap(
    gst_heat,
    color = colorRampPalette(c('deepskyblue3','grey20','firebrick1'))(11),
    #color = greenred(11),
    breaks = c(-1,-0.8,-0.6,-0.3,-0.2,0,0.2,0.3,0.6,0.8,1),
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
}



flag <- 'CXCL16'

cancer_gs_ssgsea <- list(BRCA_gs_ssgsea,THCA_gs_ssgsea,LUNG_gs_ssgsea,OV_gs_ssgsea,PRAD_gs_ssgsea,UCEC_gs_ssgsea)
names(cancer_gs_ssgsea) <-c('BRCA_gs_ssgsea','THCA_gs_ssgsea','LUNG_gs_ssgsea','OV_gs_ssgsea','PRAD_gs_ssgsea','UCEC_gs_ssgsea')

cancer_gs <- list(BRCA_gs,THCA_gs,LUNG_gs,OV_gs,PRAD_gs,UCEC_gs)
names(cancer_gs) <-c('BRCA_gs','THCA_gs','LUNG_gs','OV_gs','PRAD_gs','UCEC_gs')



setwd('~/TAM/genescoring/out/plot/CXCL16')

for (i in 1:length(cancer)){
    png(filename=paste0(names(cancer_gs_ssgsea)[[i]],'_',flag,'.png'),width=850,height=750,units = "px", bg="white" , res = 120,)
    Heat_gs(cancer_gs_ssgsea[[i]],cancer_gs[[i]],flag)
    dev.off() 
    
}



## scatter

setwd('~/TCGA_RNA_data')
BRCA2 <- as.data.frame(t(read.csv('BRCA_RNAseq_TCGA.csv',header = T, row.names = 1)))
LUNG2 <- as.data.frame(t(read.csv('LUNG_RNAseq_TCGA.csv',header = T, row.names = 1)))
PRAD2 <- as.data.frame(t(read.csv('PRAD_RNAseq_TCGA.csv',header = T, row.names = 1)))
THCA2 <- as.data.frame(t(read.csv('THCA_RNAseq_TCGA.csv',header = T, row.names = 1)))
UCEC2 <- as.data.frame(t(read.csv('UCEC_RNAseq_TCGA.csv',header = T, row.names = 1)))
OV2 <- as.data.frame(t(read.csv('OV_RNAseq_TCGA.csv',header = T, row.names = 1)))

BRCA3 <- select(BRCA2, CD163, CXCL16, CD47, CSNK2B)
LUNG3 <- select(LUNG2, CD163, CXCL16, CD47, CSNK2B)
PRAD3 <- select(PRAD2, CD163, CXCL16, CD47, CSNK2B)
THCA3 <- select(THCA2, CD163, CXCL16, CD47, CSNK2B)
UCEC3 <- select(UCEC2, CD163, CXCL16, CD47, CSNK2B)
OV3 <- select(OV2, CD163, CXCL16, CD47, CSNK2B)


BRCA3$type <- 'BRCA'
LUNG3$type <- 'LUNG'
PRAD3$type <- 'PRAD'
THCA3$type <- 'THCA'
UCEC3$type <- 'UCEC'
OV3$type <- 'OV'


scatter <- rbind(BRCA3, LUNG3, PRAD3, THCA3, UCEC3, OV3)

ggscatter(scatter, x = "CXCL16", y = "CD163", size = 0.5, 
          rug = TRUE,                                
          color = "type", palette = "jco") +
  stat_cor(aes(color = type), method = "spearman")


                  
ggscatter(scatter, x = "CXCL16", y = c("CD163", "CD47"), size = 0.3,
          combine = TRUE, ylab = "Expression",
          color = "type", palette = "jco",
          add = "reg.line", conf.int = TRUE) +
  stat_cor(aes(color = type), method = "spearman", label.x = 14)



ggscatter(scatter, x = "CXCL16", y = "CD47", size = 0.3,
          color = "type", palette = "npg",
          facet.by = "type", #scales = "free_x",
          add = "reg.line", conf.int = TRUE) +
  stat_cor(aes(color = type), method = "spearman", label.y = 2)



library(cowplot) 
# Main plot
pmain <- ggplot(scatter, aes(x = CXCL16, y = CD47, color = type))+
  geom_point()+
  ggpubr::color_palette("jco") +
  stat_cor(aes(color = type), method = "spearman")
# Marginal densities along x axis
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = scatter, aes(x = CXCL16, fill = type),
               alpha = 0.4, size = 0.2)+
  ggpubr::fill_palette("jco")
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = scatter, aes(x = CD47, fill = type),
               alpha = 0.5, size = 0.2)+
  coord_flip()+
  ggpubr::fill_palette("jco")
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)


# Scatter plot colored by groups ("Species")
sp <- ggscatter(scatter, x = "CXCL16", y = "CD47",
                color = "type", palette = "jco",
                size = 3, alpha = 0.6)+
  border()  +
  stat_cor(aes(color = type), method = "spearman")                                         
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(scatter, "CXCL16", fill = "type",
                   palette = "jco")
yplot <- ggdensity(scatter, "CD47", fill = "type", 
                   palette = "jco")+
  rotate()
# Cleaning the plots
sp <- sp + rremove("legend")
yplot <- yplot + clean_theme() + rremove("legend")
xplot <- xplot + clean_theme() + rremove("legend")
# Arranging the plot using cowplot
library(cowplot)
plot_grid(xplot, NULL, sp, yplot, ncol = 2, align = "hv", 
          rel_widths = c(2, 1), rel_heights = c(1, 2))




##############################
setwd('~/조선옥_교수님')

BRCA_clinic <- read.csv('BRCA_clinic.csv', sep = ',', header = T ,row.names = 1)

BRCA_clinic1 <- dplyr::select(BRCA_clinic, grep('.01$', names(BRCA_clinic)))

names(BRCA_clinic1) <- gsub('.01','',names(BRCA_clinic1))

BRCA <- dplyr::filter(BRCA1, BRCA1$X=='CD163' | BRCA1$X=='CD47' | BRCA1$X=='CXCL16')
rownames(BRCA) <- c('CD47','CXCL16','CD163')
BRCA <- dplyr::select(BRCA, !X)

BRCA <- dplyr::select(BRCA, grep('.01$', names(BRCA)))

names(BRCA) <- gsub('.01','',names(BRCA))



BRCA %<>% mutate_if(is.double,as.character)
t <- dplyr::bind_rows(t, TNBC)
write.csv(t,'test.csv')

TNBC <- read.csv('TNBC.csv', header = T, row.names = 1, sep = ',')

t1 <- read.csv('test.csv', sep = ',', header = T, row.names = 1)



##
T1 <- t(t1)
T1 <- as.data.frame(T1)
T1$sample <- rownames(T1)

##CD163
T1$CD163 <- as.numeric(T1$CD163)
CD163_m <- T1 %>% arrange(T1[names(T1)=='CD163'])
for (i in 1:nrow(CD163_m)){
CD163_m$pathologic_stage[[i]] <- ifelse(
  CD163_m$pathologic_stage[[i]] == 'Stage IA' || 
    CD163_m$pathologic_stage[[i]] == 'Stage I' ||
    CD163_m$pathologic_stage[[i]] == 'Stage IB' ||
    CD163_m$pathologic_stage[[i]] == 'Stage IIA' ||
    CD163_m$pathologic_stage[[i]] == 'Stage IIB' ||
    CD163_m$pathologic_stage[[i]] == 'Stage II' , 'stage A' , 
  ifelse(CD163_m$pathologic_stage[[i]] == 'Stage IIIA' ||
           CD163_m$pathologic_stage[[i]] == 'Stage III' ||
           CD163_m$pathologic_stage[[i]] == 'Stage IIIB' ||
           CD163_m$pathologic_stage[[i]] == 'Stage IV' ||
           CD163_m$pathologic_stage[[i]] == 'Stage IVA' ||
           CD163_m$pathologic_stage[[i]] == 'Stage IVB' , 'stage B' , NA))
}
CD163_m <- dplyr::filter(CD163_m, !is.na(CD163))
sample <- CD163_m$sample
CD163_m <- as.data.frame(t(CD163_m))
names(CD163_m) <- sample
CD163_m <- CD163_m[1:11,]

age <- as.data.frame(t(CD163_m))
age <- age$age_at_initial_pathologic_diagnosis
age <- as.data.frame(age)
age$age <- as.numeric(age$age)
median(na.omit(age$age))

write.csv(CD163_m,'CD163_anno.csv')




write.csv(CD163_m, 'CD163_anno.csv')

##CXCL16

T1$CXLC16 <- as.numeric(T1$CXCL16)
CXCL16_m <- T1 %>% arrange(T1[names(T1)=='CXCL16'])
for (i in 1:nrow(CXCL16_m)){
  CXCL16_m$pathologic_stage[[i]] <- ifelse(
    CXCL16_m$pathologic_stage[[i]] == 'Stage IA' || 
      CXCL16_m$pathologic_stage[[i]] == 'Stage I' ||
      CXCL16_m$pathologic_stage[[i]] == 'Stage IB' ||
      CXCL16_m$pathologic_stage[[i]] == 'Stage IIA' ||
      CXCL16_m$pathologic_stage[[i]] == 'Stage IIB' ||
      CXCL16_m$pathologic_stage[[i]] == 'Stage II' , 'stage A' , 
    ifelse(CXCL16_m$pathologic_stage[[i]] == 'Stage IIIA' ||
             CXCL16_m$pathologic_stage[[i]] == 'Stage III' ||
             CXCL16_m$pathologic_stage[[i]] == 'Stage IIIB' ||
             CXCL16_m$pathologic_stage[[i]] == 'Stage IV' ||
             CXCL16_m$pathologic_stage[[i]] == 'Stage IVA' ||
             CXCL16_m$pathologic_stage[[i]] == 'Stage IVB' , 'stage B' , NA))
}
CXCL16_m <- dplyr::filter(CXCL16_m, !is.na(CXCL16))
sample <- CXCL16_m$sample
CXCL16_m <- as.data.frame(t(CXCL16_m))
names(CXCL16_m) <- sample
CXCL16_m <- CXCL16_m[1:11,]

write.csv(CXCL16_m, 'CXCL16_anno.csv')
##CD47

T1$CD47 <- as.numeric(T1$CD47)
CD47_m <- T1 %>% arrange(T1[names(T1)=='CD47'])
for (i in 1:nrow(CD47_m)){
  CD47_m$pathologic_stage[[i]] <- ifelse(
    CD47_m$pathologic_stage[[i]] == 'Stage IA' || 
      CD47_m$pathologic_stage[[i]] == 'Stage I' ||
      CD47_m$pathologic_stage[[i]] == 'Stage IB' ||
      CD47_m$pathologic_stage[[i]] == 'Stage IIA' ||
      CD47_m$pathologic_stage[[i]] == 'Stage IIB' ||
      CD47_m$pathologic_stage[[i]] == 'Stage II' , 'stage A' , 
    ifelse(CD47_m$pathologic_stage[[i]] == 'Stage IIIA' ||
             CD47_m$pathologic_stage[[i]] == 'Stage III' ||
             CD47_m$pathologic_stage[[i]] == 'Stage IIIB' ||
             CD47_m$pathologic_stage[[i]] == 'Stage IV' ||
             CD47_m$pathologic_stage[[i]] == 'Stage IVA' ||
             CD47_m$pathologic_stage[[i]] == 'Stage IVB' , 'stage B' , NA))
}
CD47_m <- dplyr::filter(CD47_m, !is.na(CD47))
sample <- CD47_m$sample
CD47_m <- as.data.frame(t(CD47_m))
names(CD47_m) <- sample
CD47_m <- CD47_m[1:11,]

write.csv(CD47_m, 'CD47_anno.csv')

ggscatter(CD163_t, x = "CD163", y = "age_at_initial_pathologic_diagnosis", size = 0.5, add = "reg.line", conf.int = TRUE)+ stat_cor( method = "spearman")
ggscatter(CD47_t, x = "CD47", y = "age_at_initial_pathologic_diagnosis", size = 0.5, add = "reg.line", conf.int = TRUE)+ stat_cor( method = "spearman")
ggscatter(CXCL16_t, x = "CXCL16", y = "age_at_initial_pathologic_diagnosis", size = 0.5, add = "reg.line", conf.int = TRUE)+ stat_cor( method = "spearman")
CD163_t <- as.data.frame(t(CD163_m))

CD163_t <- as.data.frame(t(CD163_m))
CD163_t$CD163 <- as.numeric(CD163_t$CD163)
CD163_t$age_at_initial_pathologic_diagnosis <- as.numeric(CD163_t$age_at_initial_pathologic_diagnosis)
CD163_t$CD163 <- as.numeric(CD163_t$CD163)
t.test(as.numeric(CD163_t$CD163), as.numeric(CD163_t$age_at_initial_pathologic_diagnosis) )
CD163_t$TNBC.Subtype <- ifelse(grepl(is.na(CD163_t$TNBC.Subtype),'no_tnbc','tnbc') )
length(which(CD163_t$CD163 == 'low'))

cd163_low <-subset(CD163_t, CD163_t$CD163 == 'low')
cd163_high <-subset(CD163_t, CD163_t$CD163 == 'high')

length(which(CD163_t$age_at_initial_pathologic_diagnosis == NA))


t.test(table(CD163_t$CD163, CD163_t$TNBC.Subtype))
CD163_t$CD163 <- ifelse(CD163_t$CD163 < 9.3999, 'low' , 'high')
median(CD163_t$CD163)
CD163_t$pathologic_M[CD163_t$pathologic_M == 'MX'] <- NA
fisher.test(table(CD47_t$CD47, CD47_t$age_at_initial_pathologic_diagnosis ),simulate.p.value=TRUE,workspace=2e9)$p.value

CD47_t <- as.data.frame(t(CD47_m))
CD47_t$CD47 <- as.numeric(CD47_t$CD47)
CD47_t$age_at_initial_pathologic_diagnosis <- as.numeric(CD47_t$age_at_initial_pathologic_diagnosis)
median(CD47_t$CD47)
CD47_t$CD47 <- ifelse(CD47_t$CD47 < 11.0571, 'low' , 'high')
t.test(cd47_low$CD47, cd47_high$age_at_initial_pathologic_diagnosis )
CD47_t$pathologic_M[CD47_t$pathologic_M == 'MX'] <- NA

cd47_low <-subset(CD47_t, CD47_t$CD47 == 'low')
cd47_high <-subset(CD47_t, CD47_t$CD47  == 'high' )

t.test(CD47_t$CD47, CD47_t$age_at_initial_pathologic_diagnosis)


CXCL16_t <- as.data.frame(t(CXCL16_m))
CXCL16_t$CXCL16 <- as.numeric(CXCL16_t$CXCL16)
CXCL16_t$age_at_initial_pathologic_diagnosis <- as.numeric(CXCL16_t$age_at_initial_pathologic_diagnosis)
median(CXCL16_t$CXCL16)
CXCL16_t$CXCL16 <- ifelse(CXCL16_t$CXCL16 < 10.2914, 'low' , 'high')
chisq.test(table(CXCL16_t$CXCL16, CXCL16_t$breast_carcinoma_progesterone_receptor_status))
CXCL16_t$pathologic_M[CXCL16_t$pathologic_M == 'MX'] <- NA

cxcl_low <-subset(CXCL16_t, CXCL16_t$CXCL16 == 'low')
cxcl_high <-subset(CXCL16_t, CXCL16_t$CXCL16 == 'high')
t.test(as.numeric(CXCL16_t$CXCL16), as.numeric(CXCL16_t$age_at_initial_pathologic_diagnosis) )
