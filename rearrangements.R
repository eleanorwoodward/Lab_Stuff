source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/MafFunctions.R")

data <- read.delim("C:/Users/Noah/OneDrive/Work/Meningioma/Analysis/dRanger.txt", stringsAsFactors = F)

data3 <- ReccurentMaf(data, "gene1", 2)
PlotMaf(data3, "gene1", title = "recurrent origins of rearrangements")


data3 <- ReccurentMaf(data, "gene2", 2)
PlotMaf(data3, "gene2", title = "recurrent destination of rearrangements")

genes1 <-data[,c(1,3,24)]
genes2 <- data[,c(1,6,26)]
colnames(genes2)[2:3] <- c("chr1", "gene1")


genes.all <- rbind(genes1, genes2)
genes.all.5 <- ReccurentMaf(genes.all,"gene1", 4)
PlotMaf(genes.all.5, "gene1", title = "Rearrangements starting and ending at given gene")
PlotMaf(genes.all, "chr1", title = "Rearrangments per chromosome") 


matches <- NULL

for (i in 1:nrow(data)){
  current.gene <- data$gene1[i]
  current.target <- data$gene2[i]
  friends <- FilterMaf(data[-i,], current.gene, "gene1")
  if (nrow(friends) > 0){
      idx <- friends$gene2 %in% current.target
      close.friends <- friends[idx,]
      
      
      if (nrow(close.friends) > 0){
         temp <- c(rownames(data)[i], current.gene, current.target, data$individual[i])
         matches <- rbind(matches, temp)
         
        }
      
  } 
}

matches.table <- matrix(matches, length(matches)/4, 4)
x <- data.frame(matches.table, stringsAsFactors = FALSE)
x2 <- x
## gets rid of duplicates from same sample

x <- FilterMaf(x, c("LAMA4", "NF2", "ACO1"), "X2")

for (i in 1:nrow(x)){
  current.gene <- x[, 2][i]
  current.target <- x[, 3][i]
  current.individual <- x[, 4][i]
  friends <- FilterMaf(x[-i,], current.gene, "X2")
  if (nrow(friends) > 0){
    idx <- friends[, 3] %in% current.target
    close.friends <- friends[idx,]
    
    if (nrow(close.friends) > 0){
      for (j in 1:nrow(close.friends)){
          if (close.friends[, 4][j] == current.individual){
              x[as.numeric(rownames(close.friends)[j]), 5] <- 1
          }
              
      }
          
      }
      
      
    }
    
  } 


