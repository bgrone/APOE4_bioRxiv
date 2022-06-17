setwd("/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq")

library(plotrix)

# Make genelists with logFC and P value thresholds for the genotype topTags files
myFiles = list.files(path="/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq", pattern="^topTags_5mo_genotype_*")
myFiles <- substr(myFiles, 22,)
#myFiles <- sub('.', '', myFiles)
myFiles <- gsub("^.{0,21}", "", myFiles)

for (val in myFiles) {
    
    mo5name <- paste("topTags_5mo_genotype_", val, sep = "")
    mo10name <- paste("topTags_10mo_genotype_", val, sep = "")
    mo15name <- paste("topTags_15mo_genotype_", val, sep = "")
    mo20name <- paste("topTags_20mo_genotype_", val, sep = "")
    
    mo5 <- read.csv(mo5name)
    mo10 <- read.csv(mo10name)
    mo15 <- read.csv(mo15name)
    mo20 <- read.csv (mo20name)
    
    
    # logFC threshold
    mo5$abslogfc <- abs(mo5$logFC)
    mo10$abslogfc <- abs(mo10$logFC)
    mo15$abslogfc <- abs(mo15$logFC)
    mo20$abslogfc <- abs(mo20$logFC)
    
    mo5_sig  <- mo5[(mo5$abslogfc > 0),]
    mo10_sig  <- mo10[(mo10$abslogfc > 0),]
    mo15_sig  <- mo15[(mo15$abslogfc > 0),]
    mo20_sig  <- mo20[(mo20$abslogfc > 0),]
    
    
    # P value threshold
    mo5_sig  <- as.vector(mo5_sig[mo5_sig$FDR < 0.05 , 1])
    mo10_sig  <- as.vector(mo10_sig[mo10_sig$FDR < 0.05 , 1])
    mo15_sig  <- as.vector(mo15_sig[mo15_sig$FDR < 0.05 , 1])
    mo20_sig  <- as.vector(mo20_sig[mo20_sig$FDR < 0.05 , 1])
    
    
    # select length
    mo5_genenum  <- length(mo5_sig)
    mo10_genenum  <- length(mo10_sig)
    mo15_genenum  <- length(mo15_sig)
    mo20_genenum  <- length(mo20_sig)
    
    
    #replace blanks with NA
    mo5_sig[!length(mo5_sig) > 0] <- NA
    mo10_sig [!length(mo10_sig) > 0] <- NA
    mo15_sig [!length(mo15_sig) > 0] <- NA
    mo20_sig [!length(mo20_sig) > 0] <- NA
    
    
    genes <- data.frame(genes_5mo=NA, genes_10mo=NA, genes_15mo=NA, genes_20mo=NA) #[numeric(0), ]
    
    
    #df <- data.frame(vec1 = rep(NA, max(sapply(list(intersect, unique_mo5, unique_mo20), length))))
    genes[1:length(mo5_sig), 1] <- mo5_sig
    genes[1:length(mo10_sig), 2] <- mo10_sig
    genes[1:length(mo15_sig), 3] <- mo15_sig
    genes[1:length(mo20_sig), 4] <- mo20_sig
    
    #colnames(df) <- c("intersect", "unique_mo5", "unique_mo20")
    
    filename <- paste("genotype_pval_logFC_allgenes_", val, sep = "")
    write.csv(genes, file = filename, row.names = FALSE, na = '')
}


# Grouped Bar Plot
myFiles = list.files(path="/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq", pattern="*genotype_pval_logFC_allgenes_*")

# create data.frame for results
results <- data.frame()

for(val in myFiles)
    
{
    
    genelist <- read.csv(val)
    
    genelist[genelist==""] <- NA
    
    genes_5mo <- length(na.exclude(genelist$genes_5mo))
    genes_10mo <- length(na.exclude(genelist$genes_10mo))
    genes_15mo <- length(na.exclude(genelist$genes_15mo))
    genes_20mo <- length(na.exclude(genelist$genes_20mo))
    
    results <- rbind(results, c(genes_5mo, genes_10mo, genes_15mo, genes_20mo))
    colnames(results) <- c("5mo", "10mo", "15mo", "20mo")
    
}


#transpose, as matrix
counts <- t(results)
colnames(counts) <- rownames(results)


#plot
#cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0mo2042", "#0072B2", "#D55E00", "#CC79A7")
pdf("genotype_genelist_barplot_pval05_logfc0.pdf", width = 11)

newdata <- counts
#newdata[newdata>250]<-newdata[newdata>260]-240
barplot(newdata, main="E4vE3 Number of Differentially Expressed Genes",
        space= c(0,2),
        #ylim=c(0,300),
        #axes = FALSE,
        xlab="Cell Clusters", ylab="Genes", col=c("#E69F00", "#56B4E9", "#009E73", "#0072B2"),
        #legend = rownames(counts), 
        beside=TRUE)
#axis(2,at=c(0,50,100,150,200,275),
#     labels=c(0,50,100,150,200,575))

#box()
#axis.break(2,260,style="gap")
legend(145,350, rownames(counts), fill=c("#E69F00", "#56B4E9", "#009E73", "#0072B2"))

dev.off()
