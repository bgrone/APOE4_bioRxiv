library(VennDiagram)
library(ggplot2)

setwd("/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq/CRE/KI_CRE_difference")

#read in csv files from Taubes Supplementary Table 3
KEGG_E4inE4_vs_E3inE4 <- read.csv("/Users/briangrone/Desktop/Gladstone/Najm_Data/Najm_20190821_E4inE4_vs_E3inE4_KEGG_Excitatory.csv")
#KEGG_E4inE4_vs_E3inE3 <- read.csv("/Users/briangrone/Desktop/Gladstone/Najm_Data/Najm_20190821_E4inE4_vs_E3inE3_KEGG_Excitatory.csv")


#Extract Pathway Names
KEGG_E4inE4_vs_E3inE4 <- KEGG_E4inE4_vs_E3inE4$Description
#KEGG_E4inE4_vs_E3inE3 <- KEGG_E4inE4_vs_E3inE3$Description


selected_paths <- read.csv("/Users/briangrone/Desktop/Gladstone/Mouse_snRNAseq/CRE/KI_CRE_difference/KI_KEGG_pathway_min2_p05_summary.csv")
selected_paths <- selected_paths$X[1:16]

    #make venn diagrams
        file_name <- "Human_Mouse_Venn.pdf"
        temp1 <- venn.diagram(x = list(KEGG_E4inE4_vs_E3inE4, selected_paths), 
                     category.names = c("Human Cells", "Mouse Knockin"),
                     filename = NULL,
                     #fill = myCol,
                     height = 5, width = 5, resolution = 300, 
                     output=FALSE)
        
        ggsave(temp1, file = file_name, device = "pdf")
        
        table_name <- file_name <- paste("Human_Mouse_Venn.csv")
        overlap <- intersect(KEGG_E4inE4_vs_E3inE4, selected_paths)
        write.csv(overlap, file = table_name)
  

