library(enrichR)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(doseplot)
library(RCy3)
cytoscapePing()
library(igraph)


setwd("/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq/CRE/KI_CRE_difference")

path <- "/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq/CRE/KI_CRE_difference"
KI_files <- list.files(path=path, pattern="^topTags_5mo_genotype*")

selected_paths <- read.csv("/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq/CRE/KI_CRE_difference/KI_KEGG_pathway_min2_p05_summary.csv")
selected_paths <- selected_paths$X[1:16]

clusters <- c(1:27)

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
    
    #extract gene 
    toptags_use <- rbind(toptags_up, toptags_down)
    
    final_genes <- as.vector(na.omit(toptags_use$X))

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
                        
                        
                        heatmap2_name <- paste("KI_E4vE3_KEGG_heatmap2_5mo_nopval_", val, ".pdf", sep = "")
                        
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



