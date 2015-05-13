library(reshape2)
library(ggplot2)
library(xlsx)

#Load the data
p_vals = read.csv("H3K4me1.chr21.imputed.pval.signal.average.500bp.mx", row.names = 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, strip.white=TRUE)
cell_groups <- read.xlsx("Metadata.xlsx",1)
#Get just the group names
cell_groups <- cell_groups[,3]

pca <- prcomp(p_vals)
melted <- cbind(cell_groups, melt(pca$rotation[,1:9]))

melted <- cbind(cell_groups, melt(pca$rotation[,1:9]))
barplot <- ggplot(data=melted) +
   geom_bar(aes(x=Var1, y=value, fill=cell_groups), stat="identity") +
   facet_wrap(~Var2)
barplot