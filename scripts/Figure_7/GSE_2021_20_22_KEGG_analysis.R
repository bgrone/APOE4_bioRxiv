library(Seurat)
#library(zinbwave)
#library(edgeR)
#library(SummarizedExperiment)
library(ggplot2)
library(enrichR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ggnewscale)

setwd("/Users/briangrone/Desktop/Gladstone/Human_snRNAseq/2021_10_22")

#pathway analysis
path <- "/Users/briangrone/Desktop/Gladstone/Human_snRNAseq/2021_10_22"
human_files <- list.files(path=path, pattern="^topTags_44v33_inAD_GSE157827_2021_10_22*")

#make vector for results table
pathways <- vector()

clusters <- c(1:length(human_files))

for(val in clusters) {
    
    KI_file <- human_files[val]
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
    
    
    #pathways combining up and down genes
    #dbs <- c("GO_Biological_Process_2015") # "GO_Molecular_Function_2015",, "KEGG_2019_Mouse", "GO_Cellular_Component_2015")
    
    #pathways <- enrichr(final_genes, dbs)
    #path_name <- paste(celltype, "_20mo_knockin_minus_SynCRE_logFC0_1_pathways.csv", sep = "")
    #write.csv(pathways, path_name)
    
    genes_check <- length(final_genes)
    
    if (genes_check > 0) {
        
        #convert gene symbol to Entrez ID
        require(org.Hs.eg.db)
        ncbiID = mapIds(org.Hs.eg.db,
                        keys=final_genes,
                        column="ENTREZID",
                        keytype="SYMBOL",
                        multiVals="first")
        
        ncbi_check <- length(ncbiID)
        
        if (ncbi_check > 0) {
            ego <- enrichKEGG(
                ncbiID,
                organism = "hsa",
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
                        ego <- setReadable(ego, org.Hs.eg.db, keyType = "ENTREZID" ) # refer to man for additional parameters
                        
                        ## To color genes by log fold changes, we need to extract the log fold changes from our results table creating a named vector
                        foldchanges <- toptags_use$logFC
                        names(foldchanges) <- toptags_use$X
                        
                        #
                        celltype <- substr(KI_file, 41,nchar(KI_file)-4)
                        
                        #dotplot
                        dotplot(ego, showCategory=20, font.size=10)
                        
                        dotplot_name <- paste("44v33_inAD_GSE157827_2021_10_22_KEGG_dotplot_", celltype, ".pdf", sep = "")
                        ggsave(dotplot_name, height = 7, width = 7)
                        
                        
                        #cnetplot
                        cnetplot(ego,
                                 foldChange=foldchanges)
                        cnetplot_name <- paste("44v33_inAD_GSE157827_2021_10_22_KEGG_cnetplot_logFC_", celltype, ".pdf", sep = "")
                        ggsave(cnetplot_name, height = 9, width = 9)
                        
                        
                        ego_human <- ego@result
                        ego_human <- ego_human[ego_human$p.adjust < 0.05, ]
                        
                        csv_name <- paste("44v33_inAD_GSE157827_2021_10_22_KEGG_", celltype, ".csv", sep = "")
                        write.csv(ego_human, csv_name)
                        
                        
                        
                    }
                }
            }
        }
        
    }
}

