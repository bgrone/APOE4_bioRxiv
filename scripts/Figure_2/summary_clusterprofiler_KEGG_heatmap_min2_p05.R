library(enrichR)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(ggplot2)
library(plyr)
library(gplots)


setwd("/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq/CRE/KI_CRE_difference")

path <- "/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq/CRE/KI_CRE_difference"

ages <- c("5", "10", "15", "20")

for(age in ages) {
    
    file_pattern <- paste("^topTags_", age, "mo_genotype*", sep = "")
    
    KI_files <- list.files(path=path, pattern=file_pattern)
    #KI_files <- KI_files[1:8]
    
    #make vector for results table
    pathways <- matrix(1:500,)
    allpathways <- vector()
    
    clusters <- c(1:length(KI_files))
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
                            colnames <- cbind(colnames, paste0("cluster_", val))
                            
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
    
    
    #make heatmap
    
    file_name <- paste("summary_clusterprofiler_KI_heatmap_min2_p05_", age, "mo.pdf", sep = "")
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
    
    csv_name <- paste("summary_clusterprofiler_KI_heatmap_min2_p05_", age, "mo.csv", sep = "")
    write.csv(result_table, file = csv_name)
    
    
    
}
