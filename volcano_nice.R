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
volcano_nice <- function(df,hAss,FCIndex,pValIndex,IDIndex,vAss = NULL){
  
  df <- df[complete.cases(df),]
  names(df)[1] <- "X1"
  
  hAss <- -log(hAss)
  
  names(df)[FCIndex] <- "logFC"
  names(df)[pValIndex] <- "adj.P.Val"
  
  xlimAbs <- ceiling(max(abs(df[,FCIndex])))
  ylimAbs <- ceiling(max(abs(-log(df[,pValIndex]))))
  
  if (is.null(vAss))
  {
    vAss <- xlimAbs/10
  }
  
  xneg <- function(x) abs(hAss-1+x/(x+vAss))
  xpos <- function(x) abs(hAss-1+x/(x-vAss))
  
  test <- function(x,y,vAss){
    if (x < -vAss) {
      if (xneg(x) < -log(y)) {
        return("1")
      }
      else {
        return("0")
      }
    }
    else {
      if (x > vAss){
        if (xpos(x) < -log(y)) {
          return("1")
        }
        else {
          return("0")
        }
      }
      else {
        return("0")
      }
    }
  }
  
  df$couleur <- "0"
  df$couleur <- apply(df, 1, FUN = function(x) test(as.numeric(x[FCIndex]),as.numeric(x[pValIndex]), vAss))
  df$condLabel <- df[,IDIndex]
  df[df$couleur == "0","condLabel"] <- ""
  
  a <- ggplot(df, aes(x = logFC, y = -log(adj.P.Val), color = couleur)) + 
    geom_point(alpha = 0.5) + 
    geom_text(aes(label = condLabel, vjust = -1.5, hjust = 0.15)) +
    stat_function(fun = xneg, xlim = c(-xlimAbs,-vAss), color = "black", alpha = 0.7) + 
    ylim(c(0,ylimAbs)) + xlim(c(-xlimAbs,xlimAbs)) + 
    stat_function(fun = xpos, xlim = c(vAss,xlimAbs), color = "black", alpha = 0.7) + 
    scale_colour_manual(values = c("grey30","royalblue3")) + 
    theme_minimal() +
    theme(legend.position = "none") 
    
  return(a)
}

volcano_genesets <- function(df, gtt, consensus, outpath, mt = NULL, get_df = 0) {
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
  
  for (term in terms) 
  {
    subDf <- df[df$term == term,]
    
    if (length(subDf[,1]) > 4) 
    {
      print(term)
      print(length(subDf[,1]))
      #normalized <- (subDf$AveExpr-min(subDf$AveExpr))/(max(subDf$AveExpr)-min(subDf$AveExpr))
      a <- volcano_nice(subDf,0.05,0.5,7,6)
      print("yo")
      ggsave(paste(term,".pdf"), a, device = pdf)
    }
  }
  if (get_df == 1){
    return(df)
  }
}
