library(enrichR)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(ggplot2)
library(plyr)
library(gplots)


setwd("/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq/CRE/KI_CRE_difference")

path <- "/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq/CRE/KI_CRE_difference"

#5mo
#KI_files <- list.files(path=path, pattern="^topTags_5mo_genotype*")
#KI_files <- KI_files[c(1,2,3,6,9,10,13,17,19,21,25)]

#15mo
KI_files <- list.files(path=path, pattern="^topTags_15mo_genotype*")
KI_files <- KI_files[c(1,2,3,6,10,12,17,19,21,26,27)]


#5mo
#CRE_files <- list.files(path=path, pattern="^topTags_SynCREyoung_5mo_*")
#CRE_files <- CRE_files[c(1,2,5,11,6,8,10,13,16,19,24)]

#15mo Syn-CRE  names
CRE_files <- list.files(path=path, pattern="^topTags_Syn_CRE_E4vE3_*")
CRE_files <- CRE_files[c(1,2,3,6,10,11,14,15,17,21,23)]

#make vector for results table
pathways <- matrix(1:500,)
allpathways <- vector()

clusters <- c(1:length(CRE_files))
colnames <- vector()


for(val in clusters) {
    
    
    KI_file <- KI_files[val]
    toptags <- read.csv(KI_file)
    
    allgenes <- toptags$X
    
    # P value threshold
    toptags <- toptags[toptags$FDR < 0.05, ]
    
    #absolute value logFC threshold
    toptags_up <- toptags[toptags$logFC > 0.1, ]
    toptags_down <- toptags[toptags$logFC < -0.1, ]
    
    #extract gene names
    toptags_upgenes <- as.vector(na.omit(toptags_up$X))
    toptags_downgenes <- as.vector(na.omit(toptags_down$X))
    
    #extract genes from knockin dataset
    toptags_use <- rbind(toptags_up, toptags_down)
    
    #Syn-CRE Gene names
    CRE_file <- CRE_files[val]
    CRE_toptags <- read.csv(CRE_file)
    
    # P value threshold
    CRE_toptags <- CRE_toptags[CRE_toptags$FDR < 0.05, ]
    
    #logFC threshold
    CRE_toptags_up <- CRE_toptags[CRE_toptags$logFC > 0.1, ]
    CRE_toptags_down <- CRE_toptags[CRE_toptags$logFC < -0.1, ]
    
    #extract gene names
    CRE_toptags_upgenes <- as.vector(na.omit(CRE_toptags_up$X))
    CRE_toptags_downgenes <- as.vector(na.omit(CRE_toptags_down$X))
    
    final_upgenes <- setdiff(toptags_upgenes, CRE_toptags_upgenes)
    final_downgenes <- setdiff(toptags_downgenes, CRE_toptags_downgenes)
    final_genes <- c(final_upgenes, final_downgenes)
    
    #young_syn celltypes
    #celltypes <- c("DG1", "DG2", "CA1", "CA2_3", "Subiculum", "SST_Interneurons", "Cajal_Retzius", "Astrocytes", "OPCs", "Oligodendrocytes", "Endothelial")
    
    #15mo celltype name
    celltypes <- c("DG1", "DG2", "CA1", "CA2_3", "SST_Interneurons", "RELN_Interneurons", "Astrocytes", "OPCs", "Oligodendrocytes", "Fibroblast_like", "Choroid")
    
    celltype <- celltypes[val]
    
    genes_check <- length(final_genes)
    
    if (genes_check > 0) {
        
        #convert gene symbol to Entrez ID
        require(org.Mm.eg.db)
        ncbiID = mapIds(org.Mm.eg.db,
                        keys=final_genes,
                        column="ENTREZID",
                        keytype="SYMBOL",
                        multiVals="first")
        
        ncbi_check <- length(ncbiID)
        
        if (ncbi_check > 0) {
            ego <- enrichKEGG(
                ncbiID,
                organism = "mmu",
                minGSSize = 10,
                pvalueCutoff = 0.05,
                keyType = "ncbi-geneid" #,
                #universe = allgenes
            )
            
            
            if (!is.null(ego)) {
                
                ego <- gsfilter(ego, by='Count', min = 4)
                
                if (!is.null(ego)) {
                    
                    checkresults <- min(ego@result$p.adjust)
                    
                    if (checkresults < 0.05) {
                        
                        result  <-ego@result
                        result <- result[result$p.adjust < 0.05, ]
                        
                        pathlist <- result$Description
                        length(pathlist) <- nrow(pathways)
                        
                        #add pathways for summary table
                        pathways <- cbind(pathways, pathlist)
                        
                        colnames <- cbind(colnames, paste0("cluster_", val))
                        
                        #add pathways for total pathway table
                        allpathways  <- c(allpathways, pathlist)
                        
                    }
                }
            }
        }
        
    }
}


#make csv pathway summary table
pathways <- pathways[,-1]
colnames(pathways) <- colnames
result_table <- table(allpathways)
result_table <- result_table[result_table > 1]

result_table <- as.data.frame(result_table)
result_table <- result_table[order(result_table$Freq,  decreasing = TRUE),]
total_pathways <- result_table$allpathways

colnames(result_table) <- c("Description", "Freq")

#add p.adjust value for pathway in each cell clusters
colnames <- vector()

for(val in clusters) {
    
    KI_file <- KI_files[val]
    toptags <- read.csv(KI_file)
    
    allgenes <- toptags$X
    
    # P value threshold
    toptags <- toptags[toptags$FDR < 0.05, ]
    
    #absolute value logFC threshold
    toptags_up <- toptags[toptags$logFC > 0.1, ]
    toptags_down <- toptags[toptags$logFC < -0.1, ]
    
    #extract gene names
    toptags_upgenes <- as.vector(na.omit(toptags_up$X))
    toptags_downgenes <- as.vector(na.omit(toptags_down$X))
    
    final_genes <- c(toptags_upgenes, toptags_downgenes)
    
    genes_check <- length(final_genes)
    
    if (genes_check > 0) {
        
        #convert gene symbol to Entrez ID
        require(org.Mm.eg.db)
        ncbiID = mapIds(org.Mm.eg.db,
                        keys=final_genes,
                        column="ENTREZID",
                        keytype="SYMBOL",
                        multiVals="first")
        
        ncbi_check <- length(ncbiID)
        
        if (ncbi_check > 0) {
            ego <- enrichKEGG(
                ncbiID,
                organism = "mmu",
                minGSSize = 10,
                #pvalueCutoff = 0.05,
                keyType = "ncbi-geneid" #,
                #universe = allgenes
            )
            
            
            if (!is.null(ego)) {
                
                ego <- gsfilter(ego, by='Count', min = 4)
                
                if (!is.null(ego)) {
                    
                    checkresults <- min(ego@result$p.adjust)
                    
                    if (checkresults < 0.05) {
                        
                        pathlist <- ego@result$Description
                        length(pathlist) <- nrow(pathways)
                        
                        result <- ego@result
                        result <- result[result$p.adjust < 0.05, ]
                        
                        result_table <- join(result_table, result, by = "Description")
                        result_table <- subset(result_table, select=-c(ID, GeneRatio, BgRatio, pvalue, qvalue, geneID, Count))
                        
                        
                        #add colnames for summary table
                        colnames <- cbind(colnames, celltypes[val])
                        
                    }
                }
            }
        }
        
    }
}


colnames(result_table) <- c("Description", "Freq", as.vector(colnames))
rownames(result_table) <- result_table$Description
result_table <- subset(result_table, select=-c(Freq, Description))
result_table[result_table > 0.05] <- NA
result_table <- as.matrix(result_table)
result_table <- result_table[rowSums(is.na(result_table)) != ncol(result_table), ]
result_table <- result_table[,colSums(is.na(result_table))<nrow(result_table)]


#make heatmap
file_name <- paste("summary_clusterprofiler_KI_minus_SynCRE_heatmap_15mo.pdf", sep = "")
pdf(file = file_name     , width = 6, height = 12)

lmat = rbind(c(0,3),c(2,1),c(0,4))
lwid = c(1.5,4)
lhei = c(1.5,4,1)

heatmap.2(result_table, dendrogram = "none", trace="none", density.info = "none",
          Rowv=FALSE, Colv = FALSE,
          margin = c(4,18),
          lmat = lmat,
          lwid = lwid,
          lhei = lhei,
          cexRow = 0.8,
          cexCol = 1,
          srtCol = 45,
          na.color = "grey",
          key.title = NA, key.xlab = "P-value", key.ylab = NA
)

dev.off()

#make csv

csv_name <- paste("summary_clusterprofiler_KI_minus_SynCRE_heatmap_15mo.csv", sep = "")
write.csv(result_table, file = csv_name)




