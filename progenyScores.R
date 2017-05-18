#This is a function designed to compute progeny scores from log2foldchanges.
#Copyright (C) 2017  Aurelien Dugourd

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

progenyScores <- function(df, cm, dfIndex = 1, FCIndex = 3, cmIndex = 1) {
  names(df)[dfIndex] <- "X1"
  names(cm)[cmIndex] <- "X1"
  
  merged <- merge(df[,c(dfIndex,FCIndex)],cm)
  
  for (pathway in names(cm[,-cmIndex]))
  {
    merged[,pathway] <- merged[,2]*merged[,pathway]
  }
  
  progeny_scores <- colSums(merged[,c(3:length(merged[1,]))])
  names(progeny_scores) <- names(merged[,c(3:length(merged[1,]))])
  
  return(progeny_scores)
}



