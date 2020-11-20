library(dplyr)
library(ggplot2)
library(tidyr)
library(immunedeconv)
library(tibble)
library(knitr)
options(max.print = 99999)

setwd('~/Pink1')

df <- read.csv("pink1_norm(human_gene).csv", sep = ',', header = T)
nrow(df)
df <- na.omit(df)
df <- df[!(df$HGNC_symbol == "" ), ]

df <- arrange(df, HGNC_symbol,desc(P2),desc(P4),desc(P6),desc(P8),desc(P9),desc(W5),desc(W6),desc(W7),desc(W8),desc(W9))

df = df[-which(duplicated(df$HGNC_symbol)),] 

gene <- df$HGNC_symbol
df <- df[,-1]
row.names(df) <- gene

nrow(df)

write.csv(df, "Pink1KO_HGNC_CIBERSORT.csv")

#quantiseq

df_quantiseq = deconvolute(df, "quantiseq", tumor = F)

df_quantiseq %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_brewer(palette="Paired") +
  scale_x_discrete(limits = rev(levels(df_quantiseq)))

#mcpcounter

df_mcp_counter = deconvolute(df, "mcp_counter")

df_mcp_counter %>%
  gather(sample, score, -cell_type) %>%
  ggplot(aes(x=sample, y=score, color=cell_type)) +
  geom_point(size=4) +
  facet_wrap(~cell_type, scales="free_x", ncol=3) +
  scale_color_brewer(palette="Paired", guide=FALSE) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#knitr::kable(dataset_racle$expr_mat[1:5, ])
#knitr::kable(dataset_racle$ref[1:5, ])


#deconvolution_methods
