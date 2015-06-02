library(reshape2)
library(ggplot2)
library(xlsx)

#Pathways for files to read in
path_to_files <- paste(getwd(),"/Course_data/Epigenome_imputated_data_chrom21",sep="")
all_files <- list.files(path = path_to_files)
path_to_meta_data <- paste(getwd(), "/Course_data/Metadata.xlsx",sep="")
write_directory <- paste(getwd(), "/Course_data/Output", sep="")

cell_groups <- read.xlsx(path_to_meta_data,1)
#Get just the group names
cell_groups <- cell_groups[,3]

#10 for number of regions, 500 for size of bins, add 9 to last region (leftovers)
chr_piece <- (48129500/10)/500
chr_regions <- c(rep(1, chr_piece), rep(2, chr_piece), rep(3, chr_piece), rep(4, chr_piece), rep(5, chr_piece),
                 rep(6, chr_piece), rep(7, chr_piece), rep(8, chr_piece), rep(9, chr_piece), rep(10, chr_piece+9))



for (file in all_files){
  #Load the data and get paths for reading and writing
  file_path <- paste(path_to_files, file ,sep = "/")
  print(paste("reading file from: ", file_path, sep = ""))
  mypath <- paste(write_directory, "/PCA_BarPlot_", file, ".pdf", sep="")
  print(paste("saving barplot as: ",mypath, sep=""))
  jpeg(file=mypath)
  p_vals = read.csv(file_path, row.names = 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, strip.white=TRUE)
  new_p_vals <- p_vals[,c(1:16,18:51,53:61,63,65,66:77,79:83,87,88,90:93,98:110)]
  
  #pca
  #pca <- prcomp((p_vals))
  pca <- prcomp(t(new_p_vals))
  jpeg(paste(write_directory, "/elbow-pc-", file, ".jpg", sep=""))
  plot(pca, type='l') # figure out how to get elbow of this
  
  print(paste("Writing bins for ", file, sep=""))
  binsPath <- paste(write_directory, "/PCBins/", file, ".txt", sep="")
  write(file, file=binsPath)
  for (i in 1:5) {
    write(paste("PC", i, ":", sep=""), file=binsPath, append=TRUE)
    sorted <- sort(pca$rotation[,i], decreasing=TRUE)
    for (j in 1:length(sorted)) {
      write(paste(names(sorted)[j], sorted[j], sep=","), file=binsPath, append=TRUE)
    }
  }
  print ("Done.")
  
  #pca summary
  #x <- summary(pca)
  #vars <- x$sdev^2
  #vars <- vars/sum(vars)
  #write.table(cbind("Standard deviation" = x$sdev, "Proportion of Variance" = vars, "Cumulative Proportion" = cumsum(vars)), file=paste(write_directory, "/new-pc-", file, "-summary.txt", sep=""), sep="\t")
              
  #melted <- cbind(cell_groups, melt(pca$rotation[,1:9]))
  melted <- cbind(chr_regions, melt(pca$rotation[,1:9]))

#   #bar plot and save
#   ggplot(data=melted) +
#     geom_bar(aes(x=Var1, y=value, fill=chr_regions), stat="identity") +
#     facet_wrap(~Var2)
#   ggsave("~/Desktop/UCLA_Bioinformatics/Classes/2014/Spring/Statistical_Methods_in_Computational_Biology/Final_Project/Course_data/Output/Test3.pdf")
#   #ggsave(mypath)

  
  #pca plots
#  scores <- data.frame(cell_groups, pca$x[,1:3])
#  qplot(x=PC1, y=PC2, data=scores, colour=factor(cell_groups)) +
#    theme(legend.position="none")
#  ggsave(paste(write_directory,"/",file,"PC_12.pdf",sep=""))
#  qplot(x=PC1, y=PC3, data=scores, colour=factor(cell_groups)) +
#    theme(legend.position="none")
#  ggsave(paste(write_directory,"/",file,"PC_13.pdf",sep=""))
#  qplot(x=PC2, y=PC3, data=scores, colour=factor(cell_groups)) +
#    theme(legend.position="none")
#  ggsave(paste(write_directory,"/",file,"PC_23.pdf",sep=""))
#  print(paste("finished ", file, sep=""))

  #distance between PCs
#  d12 = dist(rbind(scores$PC1, scores$PC2))
#  print(paste("d12=", d12))
#  d23 = dist(rbind(scores$PC2, scores$PC3))
#  print(paste("d23=", d23))
#  d13 = dist(rbind(scores$PC1, scores$PC3))
#  print(paste("d13=", d13))
}