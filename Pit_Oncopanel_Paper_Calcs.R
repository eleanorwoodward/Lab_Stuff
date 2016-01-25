## Cals for Oncopanel Analysis
source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/MafFunctions.R")

## Removes trailing whitespace
whiteout <- function(string){
  x <- strsplit(string, "  ")
  x <- x[[1]]
  x <- strsplit(x, " ")
  x <- x[[1]]
  x
}

## Gets CSV data, formats and cleans it
onc.data <- read.csv('C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Data/Oncopanel_Data_For_R.csv', stringsAsFactors = FALSE, comment.char = "#", 
                   header = TRUE)
onc.data <- onc.data[1:286, ]
colnames(onc.data)[c(6, 7, 8, 9, 10, 16)] <- c("Subtype", "Recurrent", "Atypical", "MIB1_Index", "Dural_Invasion", "Gene")
onc.data$Gene <- sapply(data$Gene, whiteout, USE.NAMES = FALSE)

## Predefined clusters here
cluster.big <-c("TP53", "FANCA", "PRKDC", "ATM", "BRCA1", "ABL1", "BRCA2")
cluster.small <- c("NOTCH2", "CREBBP", "GLI3", "TCF3")
cluster.arid <- c("ARID1B", "ARID1A")
cluster4 <- c("RHPN2", "PIK3C2B")
cluster5 <- c("MEN1", "MLL2")
combo <- list(cluster.big, cluster.small, cluster.arid)
names(combo) <- c("Cluster.big", "cluser.small", "cluster.arid")
subtypes <- c("Null", "FSH", "Prolactin", "HGH", "ACTH")

fishers.helper <- function(data, group.column, group1, outcome.column, outcome1, group.names, outcome.names){
data1 <- FilterMaf(data, group1, group.column)  
data2 <- FilterMaf(data, group1, group.column, FALSE)  

data1.outcome <- FilterMaf(data1, outcome1, outcome.column)    
data2.outcome <- FilterMaf(data2, outcome1, outcome.column)    
vals <- c(nrow(data1.outcome), nrow(data1) - nrow(data1.outcome), nrow(data2.outcome), nrow(data2) - nrow(data2.outcome))
mtrx <- matrix(vals, 2, 2)
colnames(mtrx) <- group.names
rownames(mtrx) <- outcome.names
mtrx
}

## Loops through gene sets, and checks if enrichment for recurrence, atypicallity, or invasion
for(i in 1:length(combo)){
  temp <- fishers.helper(onc.data, "Gene", combo[[i]], "Recurrent", "0", c("Cluster", "Population"), c("De_Novo", "Recurrent"))
  print(temp)
  output <- fisher.test(temp)
  print(output)
}

## Loops through gene sets, and checks if enrichment for subtypes
for(i in 1:length(combo)){
  temp <- fishers.helper(onc.data, "Gene", combo[[i]], "Subtype", subtypes[4], c(names(combo)[i], "Population"), c(subtypes[4], "Population"))
  print(temp)
  output <- fisher.test(temp)
  print(output)
}

## Loops through subtypes, and checks if enrichment for recurrence, atypicallity, or invasion
for(i in 1:length(subtypes)){
  temp <- fishers.helper(onc.data, "Subtype", subtypes[[i]], "Recurrent", "0", c(subtypes[i], "Population"), c("De_Novo", "Recurrent"))
  print(temp)
  output <- fisher.test(temp)
  print(output)
}
