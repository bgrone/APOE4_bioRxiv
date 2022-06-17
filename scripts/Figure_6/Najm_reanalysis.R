library(Seurat)
library(zinbwave)
library(edgeR)
library(SummarizedExperiment)
library(ggplot2)
library(enrichR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)

setwd("/Users/briangrone/Desktop/Gladstone/Najm_Data")

load("/Users/briangrone/Desktop/Gladstone/Najm_Data/mixed_SeuratObj_broadCellType_20190821.rdata")

obj <- UpdateSeuratObject(obj)

p2 <- DimPlot(object = obj, reduction = "tsne", group.by = "celltype", shuffle = TRUE) # + NoLegend()
ggsave("mixed_SeuratObj_broadCellType_20190821_tSNE.pdf")

#select E3inE3 and E4inE4
obj <- subset(obj, sample == c("E3Ki_E3Mixed", "E4Ki_E4Mixed"))

#loop zinbwave for all clusters
clusters <- levels(obj@meta.data$celltype)

for(val in clusters)
{
    se <- subset(obj, celltype == c(val))
    
    #convert Seurat object to SummarizedExperiment
    se <- SummarizedExperiment(as.matrix(se@assays$RNA@counts),
                               colData = se@meta.data)
    
    samples <- as.integer(se@colData@nrows)
    one_percent_samples <- 0.01*samples
    
    #removing those genes that do not have more than 1 reads in more than 1% samples
    filter <- rowSums(assay(se)>1)>one_percent_samples
    table(filter)
    
    se <- se[filter,]
    
    obj_singlecell <- zinbwave(se, K=2, epsilon=1000, observationalWeights = TRUE)
    
    #differential expression
    weights <- assay(obj_singlecell, "weights")
    
    dge <- DGEList(assay(obj_singlecell))
    dge <- calcNormFactors(dge)
    
    #genotype vector for model matrix
    #se <- subset(se, se$sample %in% c("E3Ki_E3Mixed", "E4Ki_E4Mixed"))
    genotypes <- se$sample
    genotypes <- droplevels(genotypes)
    design <- model.matrix(~genotypes, data = colData(obj_singlecell))

    dge$weights <- weights
    dge <- estimateDisp(dge, design)
    fit <- glmFit(dge, design)
    
    lrt <- glmWeightedF(fit, coef = 2)
    topTags_obj <- topTags(lrt, n="all")
    file_name <- paste("topTags_Najm_20190821_E4inE4_vs_E3inE3_", val, ".csv", sep = "")
    write.csv(topTags_obj, file = file_name)
}

#pathway analysis
path <- "/Users/briangrone/Desktop/Gladstone/Najm_Data"
KI_files <- list.files(path=path, pattern="^topTags_Najm_20190821_E4inE4_vs_E3inE3_*")

#make vector for results table
pathways <- vector()

clusters <- c(1:length(KI_files))

for(val in clusters) {
    
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
                        celltype <- substr(KI_file, 40,nchar(KI_file)-4)
                        
                        #dotplot
                        dotplot(ego, showCategory=30)
                        
                        dotplot_name <- paste("Najm_20190821_E4inE4_vs_E3inE3_KEGG_dotplot_", celltype, ".pdf", sep = "")
                        ggsave(dotplot_name, height = 7, width = 7)
                        
                        
                        #cnetplot
                        cnetplot(ego,
                                 foldChange=foldchanges,
                                 showCategory=30)
                        cnetplot_name <- paste("Najm_20190821_E4inE4_vs_E3inE3_KEGG_cnetplot_logFC_", celltype, ".pdf", sep = "")
                        ggsave(cnetplot_name, height = 9, width = 9)
                        
                        
                        ego_human <- ego@result
                        ego_human <- ego_human[ego_human$p.adjust < 0.05, ]
                        
                        csv_name <- paste("Najm_20190821_E4inE4_vs_E3inE3_KEGG_", celltype, ".csv", sep = "")
                        write.csv(ego_human, csv_name)
                        
                        
                        
                    }
                }
            }
        }
        
    }
}

