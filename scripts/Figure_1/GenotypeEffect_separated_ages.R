library(zinbwave)
library(matrixStats)
library(magrittr)
library(ggplot2)
library(biomaRt)
library(edgeR)

#library(usethis)
library(SummarizedExperiment)

# Register BiocParallel Serial Execution
BiocParallel::register(BiocParallel::SerialParam())

library(Seurat)
library(dplyr)

myFiles = list.files(path="/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq", pattern="^Seurat_E3vE4.*\\.Rdata$")
myFiles = myFiles[1:27]

#data <- lapply(myFiles, read.table, sep="\t", header=FALSE)
#names(data) <- myFiles

for(val in myFiles)
    
{
    data = val
    load(data)
    
    #convert Seurat object to SummarizedExperiment for 5mo analysis
    summex <- SummarizedExperiment(as.matrix(se@assays$RNA@counts),
                                   colData = se@meta.data)
    
    #file showing number cells by age
    #cells <- table(summex@colData$age)
    #cluster = substr(data, 13, 14)
    #file_name <- paste("topTags_cellnumbers_", cluster, ".txt", sep = "")
    #write.table(cells, file = file_name)
    
    
    summex2 <- summex[ , summex@colData$age == "10"]
    
    samples <- as.integer(summex2@colData@nrows)
    one_percent_samples <- 0.01*samples
    
    #removing those genes that do not have more than 1 reads in more than 1% samples
    filter <- rowSums(assay(summex2)>1)>one_percent_samples
    table(filter)
    
    summex2 <- summex2[filter,]
    
    obj_singlecell <- zinbwave(summex2, K=2, epsilon=1000, observationalWeights = TRUE)
    
    
    #differential expression
    weights <- assay(obj_singlecell, "weights")
    
    dge <- DGEList(assay(obj_singlecell))
    dge <- calcNormFactors(dge)
    
    design <- model.matrix(~genotype, data = colData(obj_singlecell))
    dge$weights <- weights
    dge <- estimateDisp(dge, design)
    fit <- glmFit(dge, design)
    
    lrt <- glmWeightedF(fit, coef = 2)
    topTags_obj <- topTags(lrt, n="all")
    
    # save
    cluster = substr(data, 13, 14)
    file_name <- paste("topTags_10mo_genotype_", cluster, ".csv", sep = "")
    write.csv(topTags_obj, file = file_name)
    
    
    #file showing number cells by age
    #cells <- table(summex@colData$age)
    #cluster = substr(data, 13, 14)
    #file_name <- paste("topTags_cellnumbers_", cluster, ".txt", sep = "")
    #write.table(cells, file = file_name)
    
    summex3 <- summex[ , summex@colData$age == "14"]
    
    samples <- as.integer(summex3@colData@nrows)
    one_percent_samples <- 0.01*samples
    
    #removing those genes that do not have more than 1 reads in more than 1% samples
    filter <- rowSums(assay(summex3)>1)>one_percent_samples
    table(filter)
    
    summex3 <- summex3[filter,]
    
    obj_singlecell <- zinbwave(summex3, K=2, epsilon=1000, observationalWeights = TRUE)
    
    
    #differential expression
    weights <- assay(obj_singlecell, "weights")
    
    dge <- DGEList(assay(obj_singlecell))
    dge <- calcNormFactors(dge)
    
    design <- model.matrix(~genotype, data = colData(obj_singlecell))
    dge$weights <- weights
    dge <- estimateDisp(dge, design)
    fit <- glmFit(dge, design)
    
    lrt <- glmWeightedF(fit, coef = 2)
    topTags_obj <- topTags(lrt, n="all")
    
    # save
    cluster = substr(data, 13, 14)
    file_name <- paste("topTags_15mo_genotype_", cluster, ".csv", sep = "")
    write.csv(topTags_obj, file = file_name)
    
    rm(list=ls())
}
