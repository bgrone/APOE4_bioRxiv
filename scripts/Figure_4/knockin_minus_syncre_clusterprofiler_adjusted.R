library(ggplot2)
library(enrichR)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)


setwd("/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq/CRE/KI_CRE_difference")


path <- "/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq/CRE/KI_CRE_difference"
KI_files <- list.files(path=path, pattern="^topTags_15mo_genotype*")

#5mo
#KI_files <- KI_files[c(1,2,3,6,9,10,13,17,19,21,25)]

#15mo
KI_files <- KI_files[c(1,2,3,6,10,12,17,19,21,26,27)]

#15mo Syn-CRE names
CRE_files <- list.files(path=path, pattern="^topTags_Syn_CRE_E4vE3_*")

#5mo
#CRE_files <- list.files(path=path, pattern="^topTags_SynCREyoung_5mo_*")

#5mo
#CRE_files <- CRE_files[c(1,2,5,11,6,8,10,13,16,19,24)]

#15mo
CRE_files <- CRE_files[c(1,2,3,6,10,11,14,15,17,21,23)]

clusters <- c(1:11)

#make vector for results table
pathways <- matrix(1:500,)
allpathways <- vector()

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
    
    #enrichr
    #if (is.null(dbs)) websiteLive <- FALSE
    #if (websiteLive) head(dbs)
    
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
    
    #celltype name
    #15mo
    celltypes <- c("DG1", "DG2", "CA1", "CA2_3", "SST_Interneurons", "RELN_Interneurons", "Astrocytes", "OPCs", "Oligodendrocytes", "Fibroblast_like", "Choroid")
    
    #young_syn
    #celltypes <- c("DG1", "DG2", "CA1", "CA2_3", "Subiculum", "SST_Interneurons", "Cajal_Retzius", "Astrocytes", "OPCs", "Oligodendrocytes", "Endothelial")
    
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
                        
                        #dotplot
                        dotplot(ego, showCategory=10) +
                            scale_y_discrete(position = "right")
                        
                        dotplot_name <- paste("knockin_minus_SynCRE_15mo_KEGG_dotplot_adjusted_", celltype, ".pdf", sep = "") 
                        ggsave(dotplot_name, height = 7, width = 7)
                        
                        #cnetplot
                        #gene_names <- names(gene_IDs)
                        #gene_IDs <- unname(ncbiID, force = FALSE)
                        #ego@gene <- gene_names[match(ego@gene, gene_IDs)]
                        
                        ## To color genes by log fold changes, we need to extract the log fold changes from our results table creating a named vector
                        toptags_use <- toptags[toptags$X %in% final_genes, ]
                        foldchanges <- toptags_use$logFC
                        names(foldchanges) <- toptags_use$X
                        
                        ## convert gene ID to symbol
                        ego2 <- setReadable(ego, org.Mm.eg.db, keyType = "ENTREZID" ) # refer to man for additional parameters
                        
                        cnetplot(ego2, 
                                 foldChange = foldchanges,
                                 showCategory = 10)
                        
                        cnetplot_name <- paste("knockin_minus_SynCRE_15mo_KEGG_cnetplot_adjusted_", celltype, ".pdf", sep = "")
                        ggsave(cnetplot_name, height = 10, width = 10)
                        
                        pathlist <- ego2@result$Description
                        length(pathlist) <- nrow(pathways)
                        
                        #add pathways for summary table
                        pathways <- cbind(pathways, pathlist)
                        
                        colnames <- cbind(colnames, celltype)
                        
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
result_table <- as.data.frame(result_table)

#truncate pathways to match length of total resuls
pathways <- pathways[1:nrow(result_table), ]

pathways <- cbind(result_table, pathways)

#add TRUE/FALSE evaluation for presence of pathway in cell clusters
columns <- 3:ncol(pathways)

for (column in columns) { 
    presence <-  pathways[,1] %in% pathways[, column]
    pathways <- cbind(pathways, presence)
}

#remove pathways lists for each colum
pathways <- pathways[, -columns]

pathways <- ifelse(pathways[,3:ncol(pathways)] == TRUE, "X", "")
colnames(pathways) <- colnames
pathways <- cbind(result_table, pathways)

pathways <- pathways[order(-pathways$Freq),]
rownames(pathways) <- NULL

write.csv(pathways, "knockin_minus_SynCRE_clusterprofiler_KEGG_pathways_15mo.csv")

