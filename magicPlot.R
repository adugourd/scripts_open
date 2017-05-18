#This is a function designed to generate plots facilitated the visualisation of large datasets.
#Copyright (C) 2016  Aurelien Dugourd

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(ggplot2)
library(reshape)
library(pheatmap)
library(venneuler)

##Take two arguments : the dataframe and the output path

#Format of dataframe (example) :
#         Sample1 Sample2 Sample3 . . . . SampleN
#Protein1 654564  NA      5411    . . . . 5454
#Protein2 NA      NA      114     . . . . 741
#Protein3 4475    2121    NA      . . . . 12443
#.        .       .       .       . . . . .
#.        .       .       .       . . . . .
#.        .       .       .       . . . . .
#.        .       .       .       . . . . .
#.        .       .       .       . . . . .
#.        .       .       .       . . . . .
#ProteinM 5454    NA      987     . . . . 9987

#Basically you need columns as sample(or conditions) and rows as individuals(Proteins, genes, whatever..) that are measured in each condition

magicPlotMaker <- function(df, outpath){
  
  setwd(outpath)
  
  print("Dimension of dataframe :")
  print(dim(df))
  
  ##This part is just to generate the melted dataframe for ggplot
  df$ID <- row.names(df)
  
  melted_df <- melt(df)
  index <- c(1:length(melted_df[,1]))
  
  df <- df[,-length(df[1,])]
  
  #########################################
  ##                ggplots              ##
  #########################################
  
  a <- ggplot(melted_df, aes(x = value, fill = variable)) + geom_density(alpha = 0.3) + theme_minimal()
  b <- ggplot(melted_df, aes(x = index, y = value, group = variable, color = variable)) + geom_point() + theme_minimal()
  c <- ggplot(melted_df, aes(x = index, y = value, group = variable, color = variable)) + geom_boxplot() + theme_minimal()
  
  ggsave("densities.pdf", plot = a, device = pdf)
  ggsave("cloud_point.pdf", plot = b, device = pdf)
  ggsave("box_plot.pdf", plot = c, device = pdf)
  
  #########################################
  ##                PCA                  ##
  #########################################
  
  complete_df <- df[complete.cases(df),]
  print("Dimension of complete cases :")
  print(dim(complete_df))
  t_complete_df <- t(complete_df)
  PCA <- prcomp(t_complete_df ,center = TRUE, scale. = T) 

  pdf("boulder_plot.pdf")
  plot(PCA)
  dev.off()
  pdf("PCA.pdf")
  plot(PCA$x[,1],PCA$x[,2], pch = 19, xlab = paste("PC1 (",as.character(round(PCA$sdev[1]^2/sum(PCA$sdev^2)*100)),"%)"), ylab = paste("PC2 (",as.character(round(PCA$sdev[2]^2/sum(PCA$sdev^2)*100)),"%)"), main = "PCA of samples") 
  text(PCA$x[,1],PCA$x[,2], labels = names(PCA$x[,1]), pos = 3)
  dev.off()
  
  #########################################
  ##                heatmap              ##
  #########################################
  
  pheatmap(t_complete_df, show_colnames =  FALSE, filename = "profil_heatmap.pdf")
  pheatmap(cor(complete_df, method = "spearman"), filename = "correlation_heatmap.pdf")
}
