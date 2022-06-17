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

df[is.na(df)] = 0

#make dataframe fraction of total
df_percentages <- df / rowSums(df) * 100

df_rows <- rownames(df_percentages)

for (level in groups) {
    
    sample_list <- samples_by_group[,level]

    df_rows <- replace(df_rows, df_rows == sample_list[1], level)
    df_rows <- replace(df_rows, df_rows == sample_list[2], level)
    df_rows <- replace(df_rows, df_rows == sample_list[3], level)
    df_rows <- replace(df_rows, df_rows == sample_list[4], level)
    
}

df_rows <- sub("14mo", "15mo", df_rows)

df_rows <- sub("20mo", "mo20", df_rows)
df_rows <- sub("15mo", "mo15", df_rows)
df_rows <- sub("10mo", "mo10", df_rows)
df_rows <- sub("5mo", "mo05", df_rows)

df_percentages$group <- df_rows
df_percentages$geno <- substr(df_percentages$group, 1,2)
df_percentages$age <- substr(df_percentages$group, 4,7)

df_percentages$group <- NULL

clusters <- 1:27
for (cluster in clusters) {

cols <- c(cluster, 28, 29)

select_vals <- df_percentages[, cols]
colnames(select_vals) <- c("samples", "age", "geno")

res.aov2 <- aov(samples ~ age * geno, data = select_vals)

file_name <- paste("Cell_Count_ANOVA_", cluster, ".txt", sep = "")

#export it to wordfile
capture.output(summary(res.aov2),file=file_name)

}
