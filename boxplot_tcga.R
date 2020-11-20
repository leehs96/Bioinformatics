getwd()
setwd('../')
library(reshape)
library(ggpubr)
library(rstatix)
library(dplyr)
getwd()
setwd('../')

library(readr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(plyr)
library(RColorBrewer)
library(reshape2)

source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
          
my_data<-read.csv(url('https://data.bris.ac.uk/datasets/112g2vkxomjoo1l26vjmvnlexj/2016.08.14_AnxietyPaper_Data%20Sheet.csv'))
                                   
head(X)
df <- read.csv('SKCM_RNAseq_TCGA.csv', sep = ',', header = T, row.names = 1)
df <- read.table('pancaner.txt', sep = '\t', header = T)
#df <- read.csv('LUNG_RNAseq_TCGA.csv', sep = ',', header = T, row.names = 1)
#df <- df[names(df) == 'CD14',]
#write.csv(df,'pancancer.csv')


pri <- select(df, grep('.01$', names(df)))
meta <- select(df, grep('.06$', names(df)))


pri_CD14 <- pri[rownames(pri) == 'CD14',]
meta_CD14 <- meta[rownames(meta) == 'CD14',]

pri_CD14_1 <- melt(as.matrix(pri))
meta_CD14_1 <-melt(as.matrix(meta))

box <- rbind(pri_CD14_1, meta_CD14_1)
names(box) <- c('gene','type','value')

box$type <- ifelse(grepl(".01", box$type),"Primary tumor","Metastatic tumor")

box_CD14 <- ggviolin(box,
                      'type',
                      'value',
                      fill = 'type',
                      palette = c("#00798c","#d1495b")
                      ) + 
  xlab('Pan-cancer') +
  ylab('CD14 expression') +
  #geom_jitter(alpha = 0.15, color = "#2e4057" , position=position_jitter(0.5), shape=16) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1) +
  geom_boxplot(width = 0.4)


stat.test_CD14 <- box %>% t_test(value ~ type) %>% adjust_pvalue() %>% add_significance("p.adj")
stat.test_CD14 <- stat.test_CD14 %>% add_xy_position(x = "type")
stat.test_CD14


box_CD14 + stat_pvalue_manual(stat.test_CD14, label = "T-test, p = {p.adj}", vjust = -1, bracket.nudge.y = 2) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))


raincloud_theme = theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(angle = 45, vjust = 0.5),
  legend.title=element_text(size=16),
  legend.text=element_text(size=16),
  legend.position = "right",
  plot.title = element_text(lineheight=.8, face="bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))







g <- ggplot(box, aes(y = value, x = type, fill = type)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = value, color = type), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  #guides(fill = FALSE) +
  #guides(color = FALSE) +
  scale_color_brewer(palette = 'Dark2') +
  scale_fill_brewer(palette = "Dark2") +
  # coord_flip() +
  theme_bw() +
  raincloud_theme +  theme(legend.position="none")

g
