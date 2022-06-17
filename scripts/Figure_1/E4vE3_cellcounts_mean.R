library(dplyr)

setwd("/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq")

myFiles = list.files(path="/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq", pattern="^Seurat_E3vE4.*\\.Rdata$")
#myFiles = myFiles[8]

#data <- lapply(myFiles, read.table, sep="\t", header=FALSE)
#names(data) <- myFiles

#get full list of samples
load(myFiles[1])
allcells <- table(se@meta.data$sampleNumber$x)
allsamples <- names(allcells)

#get sample IDs for each group
se@meta.data$group <- se@meta.data$group$V1
groups <- levels(se@meta.data$group)

samples_by_group <- data.frame()

for (level in groups) {
    
    cell_subset <- subset(se, subset = group == level)
    samples <- names(table(cell_subset@meta.data$sampleNumber$x))
    
    samples_by_group <- rbind(samples_by_group, samples)
    
}

rownames(samples_by_group) <- groups
samples_by_group <- t(samples_by_group)
rownames(samples_by_group) <- NULL

#create df and name columns
df <- data.frame(matrix(ncol = 0, nrow = 32))
rownames(df) <- names <- names(allcells)

length(myFiles)

for(val in myFiles)
    
{
    data = val
    load(data)
    
    #
    cells <- table(se@meta.data$sampleNumber$x)
    samplenames <- names(cells)
    cells <- as.vector(cells)
    names(cells) <- samplenames
    #cells <- t(cells)
    
    cells <- as.data.frame(cells)

    df <- merge(df, cells, by = 0, all = TRUE)
    rownames(df) <- df$Row.names
    df$Row.names <- NULL
}


#name data.frame
colnames(df) <- 1:27

#make dataframe fraction of total
rowSums(df, na.rm = TRUE)

df_percentages <- df / rowSums(df, na.rm = TRUE) * 100

#make data frame with mean
mean_df <- data.frame()
se_df <- data.frame()

for (level in groups) {

sample_list <- samples_by_group[,level]

df_subset <- df_percentages[rownames(df_percentages) %in% sample_list, ]

df_subset[is.na(df_subset)] = 0

df_subset_mean <- colSums(df_subset)/4

mean_df <- rbind(mean_df, df_subset_mean)

df_subset_sd <- apply(df_subset, 2, sd)
df_subset_se <- df_subset_sd/sqrt(4)
se_df <- rbind(se_df, df_subset_se)

}

#finish mean matrix
colnames(mean_df) <- 1:27
rownames(mean_df) <- groups
mean_df <- as.matrix(mean_df)

#name data.frame
names <- rownames(mean_df)
names <- sub("14mo", "15mo", names)
rownames(mean_df) <- names

#reorder groups
mean_df <- mean_df[c(4,1,2,3,8,5,6,7), ]

#SE dataframe
colnames(se_df) <- 1:27
rownames(se_df) <- groups
se_df <- as.matrix(se_df)

#name data.frame
names <- rownames(se_df)
names <- sub("14mo", "15mo", names)
rownames(se_df) <- names

#reorder groups
se_df <- se_df[c(4,1,2,3,8,5,6,7), ]


#plot
#cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0mo2042", "#0072B2", "#D55E00", "#CC79A7")
pdf("E4vE3_cellcounts_mean.pdf", width = 18)

barCenters <- barplot(mean_df, 
        space= c(0,3),
        ylim = c(0,27),
        xlab="Cell Clusters", ylab="percentage of cells for each genotype and age group in each cluster", col = c("#D2FFFE", "#C1EEEE", "#A6CBCC", "#728C8C", "#E68A6B", "#D58063", "#BA7058", "#7D4C3A"),
        legend = rownames(mean_df), 
        args.legend = list(x = "topright", inset = c(0.05, 0)),
        beside=TRUE)

segments(barCenters, mean_df - se_df, barCenters,
         mean_df + se_df, lwd = 1.5)

#arrows(barCenters, mean_df - se_df, barCenters,
#       mean_df + se_df, lwd = 1.5, angle = 90,
#       code = 3, length = 0.05)

dev.off()

