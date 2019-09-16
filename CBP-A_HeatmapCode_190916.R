## RStudio merging and heatmap code for publication:

# Journal: Comparative Biochemistry and Physiology Part A: Molecular and Integrative Physiology
# Authors: Andrew J. Cline, Dr. Scott L. Hamilton, Dr. Cheryl A. Logan
# Title: Effects of multiple climate change stressors on gene expression in blue rockfish (Sebastes mystinus)

install.packages('tidyverse')
install.packages('reshape2')
install.packages('ggplot2')
install.packages('ggdendro')
install.packages('scales')

library(tidyverse) 
library(reshape2)  
library(ggplot2)   
library(ggdendro)
library(scales)  

# Merge the log2FC gene expression data with the differentially expressed contigs from each subset file:

twelve=read.csv("diffExpr.P0.05_C0.matrix.log2.centered.dat_12h",sep="\t")

twelve_pH=read.csv("matrix.counts.matrix.12h_control_vs_12h_pH.edgeR.DE_results.P0.05_C0.DE.subset",sep="\t")

twelve_pH_merged=right_join(twelve, twelve_pH, by="contig")

colnames(twelve_pH_merged)  ## to check which columns have your log2FC data (not n=4 for all...)

twelve_pH_merged_FC=twelve_pH_merged[,1:17]

twelveh_pH_merged_FC=twelve_pH_merged_FC[,2:length(twelve_pH_merged_FC)]

rownames(twelveh_pH_merged_FC)=twelve_pH_merged_FC[,1]

row.order=hclust(dist(twelveh_pH_merged_FC))$order

col.order=c(1,2,3,4,5,6,7,8) # THIS IS WHAT HOW YOU CHANGE WHICH SAMPLES ARE IN THE HEATMAP

twelveh_pH_merged_FC_new=twelveh_pH_merged_FC[row.order,col.order]

clust_twelveh_pH_merged_FC_new=melt(as.matrix(twelveh_pH_merged_FC_new))

names(clust_twelveh_pH_merged_FC_new)[c(1:2)] = c("contig","treatment")

g<-ggplot(clust_twelveh_pH_merged_FC_new, aes(x=treatment,y=contig)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low="blue", mid="black", high="yellow", midpoint=0,    limits=c(-2,2), oob=squish) + # change color scale min/max
  ylab("") +
  xlab("") +
  ggtitle("") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=10),
        axis.title = element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(angle = 20, size = 5))+
  labs(fill = "Log2 Fold Change")

g