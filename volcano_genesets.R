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

#df = result from limma differnecial expression (with column X as identifier column)
#gtt = gene to term table
#consensus = consensus table of piano analisys (with column X as term column)
#outpath = path to the directory were the figures will be stored
#mt = mapping table between the identifiers of limma analisys and gene to term table

volcano_genesets <- function(df, gtt, consensus, outpath, mt = NULL) {

gtt <- as.data.frame(apply(gtt, 2, function(x) toupper(x)))
names(gtt) <- c("gene","term")
colnames(df)[1] <- "X"
colnames(consensus)[1] <- "X"
if (!is.null(mt))
{
  mt <- apply(mt, 2, function(x) toupper(x))
  
  names(mt) <- c("X","gene")

  gtt <- merge(gtt,mt, all = FALSE)
  
  names(gtt) <- c("gene","term","X")
}
else
{
  names(gtt) <- c("X","term")
}

df$X <- toupper(df$X)
df <- merge(df,gtt, all = FALSE)

setwd(outpath)

terms <- consensus$X

df$baseMean <- log(df$baseMean)

df$normalized <- (df$baseMean-min(df$baseMean))/(max(df$baseMean)-min(df$baseMean))

for (term in terms) 
{
  subDf <- df[df$term == term,]

  if (length(subDf[,1]) > 1) 
  {
    print(term)
    print(length(subDf[,1]))
    #normalized <- (subDf$baseMean-min(subDf$baseMean))/(max(subDf$baseMean)-min(subDf$baseMean))
    subDf_005 <- subset(subDf, padj <= 0.05)
    position = (max(-log(subDf$padj))/80)+3
    #subDf$y <- (1/(sqrt(abs(subDf$log2FoldChange))))*5
    
    if (length(subDf_005[,1]) > 0)
    {
      a <- ggplot(subDf) + 
        geom_text(data = subDf_005, 
                  aes(x = log2FoldChange,y = -log(padj), label = X, vjust = -0.1, hjust = -0.35, nudge_y = 0, check_overlap = TRUE), size = 7) +
        geom_point(aes(x=log2FoldChange,y=-log(padj), fill = -log(padj)),pch = 21, size = subDf$normalized*15) + 
        #theme(legend.position = "none") +
        theme(legend.position = "none", 
              axis.title.x=element_blank(),
              axis.title.y=element_blank()) + 
        geom_vline(xintercept = 0.5, alpha = 0.3 ) + geom_vline(xintercept = -0.5, alpha = 0.3 ) +
        geom_hline(yintercept = 3, alpha = 0.3 ) +
        #annotate("text", x = -5, y = position, label = "pValue = 0.05") +
        scale_x_continuous(breaks = round(seq(-100, 100, by = 0.5),1)) + scale_fill_gradient(low="black", high="green") + theme_bw()
    }
    else
    {
      a <- ggplot(subDf) + 
        geom_point(aes(x=log2FoldChange,y=-log(padj), fill = -log(padj)),pch = 21, size = subDf$normalized*15) + 
        #theme(legend.position = "none") +
        theme(legend.position = "none", 
              axis.title.x=element_blank(),
              axis.title.y=element_blank()) + 
        geom_vline(xintercept = 0.5, alpha = 0.3 ) + geom_vline(xintercept = -0.5, alpha = 0.3 ) +
        geom_hline(yintercept = 3, alpha = 0.3 ) +
        #annotate("text", x = -5, y = position, label = "pValue = 0.05") +
        scale_x_continuous(breaks = round(seq(-100, 100, by = 0.5),1)) + scale_fill_gradient(low="black", high="green") + theme_bw()
    }
    ggsave(paste(term,".pdf"), a, device = pdf)
  }
}
}

volcano_genesets(tableTop, gene_to_term, consensus_hm, "~/Documents/Thorsten/Data_visualisation/piano_volcano_hm/", mapping_unip_to_geneName)
