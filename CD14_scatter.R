

path <- "/users/hslee/TCGA_RNA_data/oh/M"
data2 <-list.files(path = path, pattern = '*.csv', full.names = T)
############################################
for (i in seq_along(data2)){
  a <- gsub(path,"",data2[i])
  b <- gsub('/','',a)
  name <- (gsub("_RNAseq_TCGA.csv","_M_scatter",b))
  a <- (read.csv(data2[i], header = T , sep = ',', row.names = 1))
  a <- dplyr::select(a, grep('.06$', names(a)))
  a <- as.data.frame(t(a))
  assign(name, a)
}

setwd('~/TCGA_RNA_data')
normal <- t(read.csv('GTEX_TCGA_gene.csv', sep = ',', header = T, row.names = 1))
normal_ID <- read.csv('GTEx_Analysis_v8_Annotations_SampleAttributesDS.csv', sep = ',', header = T, row.names = 1 )

normal <- as.data.frame(normal)

normal_ID_CD14 <- normal_ID[grep("abnormal", normal_ID$SMPTHNTS),]
grep("Thyoird|Breast|Skin|Kidney|Lung", normal_ID_CD14$SMTS)


BRCA_scatter$type <- 'BRCA'
SKCM_M_scatter$type <- 'SKCM_meta'
SKCM_scatter$type <- 'SKCM'
THCA_scatter$type <- 'THCA'
KIRP_scatter$type <- 'KIRP'
KIRC_scatter$type <- 'KIRC'
KICH_scatter$type <- 'KICH'
LUNG_scatter$type <- 'LUNG'

scatter <- rbind(BRCA_scatter,
                 SKCM_M_scatter,
                 SKCM_scatter,
                 THCA_scatter,
                 KIRP_scatter,
                 KIRC_scatter,
                 KICH_scatter,
                 LUNG_scatter)


ggscatter(scatter, x = "CD14", y = "CD68", size = 0.5, 
          rug = TRUE,                                
          color = "type", palette = "jco") +
  stat_cor(aes(color = type), method = "spearman")



ggscatter(scatter, x = "CXCL16", y = c("CD163", "CD47"), size = 0.3,
          combine = TRUE, ylab = "Expression",
          color = "type", palette = "jco",
          add = "reg.line", conf.int = TRUE) +
  stat_cor(aes(color = type), method = "spearman", label.x = 14)



ggscatter(scatter, x = "CD14", y = "SIRPA", size = 0.3,
          color = "type", palette = "npg",
          facet.by = "type", #scales = "free_x",
          add = "reg.line", conf.int = TRUE) +
  stat_cor(aes(color = type), method = "spearman", label.y = 2)



library(cowplot) 
# Main plot
pmain <- ggplot(scatter, aes(x = CD14, y = SIRPA, color = type))+
  geom_point()+
  ggpubr::color_palette("jco") +
  stat_cor(aes(color = type), method = "spearman")
# Marginal densities along x axis
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = scatter, aes(x = CD14, fill = type),
               alpha = 0.4, size = 0.2)+
  ggpubr::fill_palette("jco")
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = scatter, aes(x = SIRPA, fill = type),
               alpha = 0.5, size = 0.2)+
  coord_flip()+
  ggpubr::fill_palette("jco")
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)


# Scatter plot colored by groups ("Species")
sp <- ggscatter(scatter, x = "CD14", y = "CD68",
                color = "type", palette = "jco",
                size = 3, alpha = 0.6)+
  border()  +
  stat_cor(aes(color = type), method = "spearman")                                         
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(scatter, "CD14", fill = "type",
                   palette = "jco")
yplot <- ggdensity(scatter, "CD68", fill = "type", 
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





