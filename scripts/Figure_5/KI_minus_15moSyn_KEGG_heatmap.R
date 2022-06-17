library(enrichR)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(doseplot)
library(RCy3)
cytoscapePing()
library(igraph)
library(gplots)


setwd("/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq/CRE/KI_CRE_difference")

selected_paths <- read.csv("/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq/CRE/KI_CRE_difference/KI_KEGG_pathway_min2_p05_summary.csv")
selected_paths <- selected_paths$X[1:16]

path <- "/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq/CRE/KI_CRE_difference"

KI_files <- list.files(path=path, pattern="^topTags_15mo_genotype*")
KI_files <- KI_files[c(1,2,3,6,10,12,17,19,21,26,27)]

#Syn-CRE Gene names
CRE_files <- list.files(path=path, pattern="^topTags_Syn_CRE_E4vE3_*")
CRE_files <- CRE_files[c(1,2,3,6,10,11,14,15,17,21,23)]

clusters <- c(1:length(CRE_files))

for(val in clusters) {
    
    #make vector for results table
    pathways <- vector()
    
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
    
    intersect_downgenes <- intersect(toptags_downgenes, CRE_toptags_downgenes)
    write.csv(intersect_downgenes, file = "intersect_KI_SYNCre_15mo_CA1_downgenes.csv")
    
    #celltype name
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
                        
                        ego <- pairwise_termsim(ego)
                        
                        ## convert gene ID to symbol
                        ego <- setReadable(ego, org.Mm.eg.db, keyType = "ENTREZID" ) # refer to man for additional parameters
                        
                        ## To color genes by log fold changes, we need to extract the log fold changes from our results table creating a named vector
                        foldchanges <- toptags_use$logFC
                        names(foldchanges) <- toptags_use$X
                        foldchanges <- foldchanges[ names(foldchanges) %in% final_genes]
                        
                        #extract results to make heatmap2 plot
                        ego_result <- ego@result
                        #ego_result <- ego_result[ego_result$p.adjust < 0.05,]
                        
                        m <- matrix(0, ncol = length(names(foldchanges)), nrow = nrow(ego_result))
                        rownames(m) <- ego_result$Description
                        colnames(m) <- names(foldchanges)
                        
                        
                        Pathways <- 1:length(rownames(m))
                        Genes <- colnames(m)
                        
                        for (pathway in Pathways) {
                            
                            for (gene in Genes) {
                                
                                if (grepl(gene, ego@result$geneID[pathway])) {
                                    
                                    m[pathway, gene]  <- unname(foldchanges[gene]) 
                                    
                                    
                                }
                                
                            }
                            
                        }
                        
                        
                        #select pathways in top shared pathways
                        m1 <- as.data.frame(m)
                        m1$pathnames <- rownames(m1)
                        
                        m1 <- m1[m1$pathnames %in% selected_paths, ]
                        
                        #remove pathnames column
                        m1 <- subset(m1, select=-c(pathnames))
                        m1 <- as.matrix(m1)
                        
                        #remove genes not detected in selected pathways
                        m1 =  m1[,colSums(m1) != 0]
                        
                        
                        heatmap2_name <- paste("KI_minus_SynCRE_E4vE3_KEGG_heatmap2_15mo_", celltype, ".pdf", sep = "")
                        
                        pdf(file = heatmap2_name, width = 48, height = 16)
                        
                        #colors = c(seq(-0.6,-.000001,length=100),seq(0,0,length=100), seq(0.000001,6,length=100))
                        
                        #my_palette <- colorRampPalette(c("blue", "white","red"))(n = 299)
                        
                        par(oma=c(5,2,2,40))
                        heatmap.2(m1, 
                                  density.info="none",  # turns off density plot inside color legend
                                  trace="none",         # turns off trace lines inside the heat map
                                  #margins =c(14,32),     # widens margins around plot
                                  #key = FALSE,
                                  #Rowv = NULL,
                                  lhei = c(1.5,13), # 
                                  #col= my_palette,       # use on color palette defined earlier
                                  #breaks=col_breaks,    # enable color transition at specified limits
                                  col = bluered,
                                  cexRow = 3,
                                  cexCol = 2,
                                  symbreaks = TRUE,
                                  dendrogram="both")            # turn off column clustering
                        
                        dev.off()
                        
                        
                        
                        
                    }
                }
            }
        }
        
    }
}



