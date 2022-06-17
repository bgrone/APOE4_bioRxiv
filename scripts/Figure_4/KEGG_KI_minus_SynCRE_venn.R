library(enrichR)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(ggplot2)
library(plyr)
library(gplots)
library(VennDiagram)
library(RColorBrewer)


setwd("/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq/CRE/KI_CRE_difference")
#Make Venn Diagram across ages


#read in csv files from E4vE3 pathways
KEGG_5mo <- read.csv("/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq/CRE/KI_CRE_difference/summary_clusterprofiler_KI_minus_SynCRE_heatmap_5mo.csv")
KEGG_10mo <- read.csv("/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq/CRE/KI_CRE_difference/summary_clusterprofiler_KI_minus_SynCRE_heatmap_10mo.csv")
KEGG_15mo <- read.csv("/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq/CRE/KI_CRE_difference/summary_clusterprofiler_KI_minus_SynCRE_heatmap_15mo.csv")


#Extract Pathway Names
KEGG_5mo <- KEGG_5mo$X
KEGG_10mo <- KEGG_10mo$X
KEGG_15mo <- KEGG_15mo$X

#make venn diagrams

file_name <- paste("KEGG_KI_minus_SynCRE_min2_p05_Venn2.pdf")
myCol <- brewer.pal(3, "Pastel1")

temp1 <- venn.diagram(x = list(Five_Month = KEGG_5mo,
                               Ten_Month = KEGG_10mo,
                               Fifteen_Month = KEGG_15mo
                               
), 
#category.names = c(genotype, "mouse_KI"),
filename = NULL,
fill = myCol,
height = 5, width = 5, resolution = 300, 
output=FALSE)


ggsave(temp1, file = file_name, device = "pdf")
#pdf(file = file_name, width = 6, height = 6)
#temp1
#dev.off()

#list table
allpaths <- unique(c(KEGG_5mo, KEGG_10mo, KEGG_15mo))

paths <- data.frame(allpaths,1)
paths$mo5 <- allpaths %in% KEGG_5mo
paths$mo10 <- allpaths %in% KEGG_10mo
paths$mo15 <- allpaths %in% KEGG_15mo

rownames(paths) <- paths$allpaths
paths <- paths[,-c(1,2)]

paths$sum <- rowSums(paths)
paths <- paths[order(paths$sum, decreasing = TRUE), ]
paths <- ifelse(paths[,1:3] == TRUE, "X", "")

write.csv(paths, "KI_minus_SynCRE_KEGG_pathway_min2_p05_summary.csv")

